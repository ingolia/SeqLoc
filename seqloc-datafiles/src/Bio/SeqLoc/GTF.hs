{-# LANGUAGE BangPatterns, MagicHash, OverloadedStrings, UnboxedTuples #-}

{-| Utilities for reading and writing GTF format gene annotations -}

module Bio.SeqLoc.GTF
       ( readGtfTranscripts, readGtfTranscriptsErr
       , transcriptToGtf
       )
       where 

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.ByteString.Internal (c2w)
import Data.Char
import Data.List
import Data.Maybe

import qualified Data.Attoparsec.Char8 as AP (isSpace_w8)
import qualified Data.Attoparsec.Zepto as ZP
import qualified Data.Conduit as C
import qualified Data.Conduit.List as C
import qualified Data.HashMap.Strict as HM
import qualified Data.Iteratee as Iter
import qualified Data.Iteratee.Char as IterChar
import qualified Data.Vector as V

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.SeqLoc.ZeptoUtils

-- | Convert a 'Transcript' to a string consisting of GTF lines. These
-- lines will contain @exon@ lines for the transcript, as well as
-- @CDS@ lines if the 'Transcript' has a 'cds'.
transcriptToGtf :: BS.ByteString -> Transcript -> BS.ByteString
transcriptToGtf src trx = BS.unlines $ exonLines ++ cdsLines
  where exonLines = splocLines src trx "exon" (location trx)
        cdsLines = maybe [] (splocLines src trx "CDS") $! cdsLocation trx
        
splocLines :: BS.ByteString -> Transcript -> BS.ByteString -> SpliceSeqLoc -> [BS.ByteString]
splocLines src trx ftype (OnSeq refname sploc) = map contigLines . Loc.toContigs $ sploc
  where contigLines contig = let (start0, end0) = Loc.bounds contig
                                 strchr = case Loc.strand contig of 
                                   Plus -> "+"
                                   Minus -> "-"
                             in unfields [ unSeqLabel refname
                                         , src
                                         , ftype
                                         , BS.pack . show . Pos.unOff . (+ 1) $ start0
                                         , BS.pack . show . Pos.unOff . (+ 1) $ end0
                                         , "0.0"
                                         , strchr
                                         , "."
                                         , attrs
                                         ]
        attrs = BS.concat [ "gene_id \"", unSeqLabel . geneId $ trx
                          , "\"; transcript_id \"", unSeqLabel . trxId $ trx
                          , "\"; " ]
        unfields = BS.intercalate (BS.singleton '\t')
                
-- | Read a GTF annotation file. The entire file is read at once,
-- because a single annotated transcript can span many lines in a GTF
-- file that are not required to occur in any specific order. The
-- transcript 'SpliceSeqLoc' transcript location is assembled from
-- @exon@ annotations and any CDS location is then produced from @CDS@
-- annotations, with an error occurring if the CDS is not a single
-- contiguous location within the transcript.
readGtfTranscripts :: FilePath -> IO [Transcript]
readGtfTranscripts = Iter.fileDriver gtfTrxsIter >=> 
                     either (ioError . userError) return . mkTranscripts
  where gtfTrxsIter = Iter.joinI . IterChar.enumLinesBS . Iter.joinI . gtflineIter $ Iter.foldl' insertGtfLine trxs0
        trxs0 = GtfTrxs HM.empty HM.empty HM.empty

readGtfTranscriptsErr :: FilePath -> IO ([Transcript], [String])
readGtfTranscriptsErr = Iter.fileDriver gtfTrxsIter >=>
                        return . mkTranscriptsErr
  where gtfTrxsIter = Iter.joinI . IterChar.enumLinesBS . Iter.joinI . gtflineIter $ Iter.foldl' insertGtfLine trxs0
        trxs0 = GtfTrxs HM.empty HM.empty HM.empty                        

mkTranscripts :: GtfTrxs -> Either String [Transcript]
mkTranscripts trxs = go [] allTrxs
  where allTrxs = HM.toList . gtfTogene $ trxs
        go curr [] = Right curr
        go curr (trxAndGene:rest) = case mkOne trxAndGene of
          Left err -> Left err
          Right t -> let next = t : curr in next `seq` go next rest
        mkOne (trxname, genename) = mkTranscript trxname exons cdses genename
          where exons = fromMaybe [] . HM.lookup trxname . gtfExonLocs $ trxs
                cdses = fromMaybe [] . HM.lookup trxname . gtfCdsLocs $ trxs

mkTranscriptsErr :: GtfTrxs -> ([Transcript], [String])
mkTranscriptsErr trxs = go ([], []) allTrxs
  where allTrxs = HM.toList . gtfTogene $ trxs
        go curr [] = curr
        go (currtrx, currerr) (trxAndGene:rest) = case mkOne trxAndGene of
          Left err -> let nexterr = err : currerr
                      in nexterr `seq` go (currtrx, nexterr) rest
          Right t -> let next = t : currtrx
                     in next `seq` go (next, currerr) rest
        mkOne (trxname, genename) = mkTranscript trxname exons cdses genename
          where exons = fromMaybe [] . HM.lookup trxname . gtfExonLocs $ trxs
                cdses = fromMaybe [] . HM.lookup trxname . gtfCdsLocs $ trxs        

mkTranscript :: BS.ByteString -> [ContigSeqLoc] -> [ContigSeqLoc] -> BS.ByteString -> Either String Transcript
mkTranscript trx exons cdses gene = moderr $ do loc <- exonLoc exons
                                                cdsloc <- cdsLoc loc cdses
                                                return $ Transcript (toSeqLabel gene) (toSeqLabel trx) loc cdsloc
  where moderr = either (Left . (("Transcript " ++ show trx ++ ": ") ++)) Right
                                                
exonLoc :: [ContigSeqLoc] -> Either String SpliceSeqLoc
exonLoc exons = do
  (seqname, rawcontigs) <- allSameName exons
  contigs <- sortclocs rawcontigs
  sploc <- maybe (badClocs contigs) Right $! SpLoc.fromContigs contigs
  return $! OnSeq seqname sploc
  where sortclocs clocs = maybe (badClocs clocs) Right . sortContigs $ clocs
        badClocs clocs = Left . unwords $ [ "bad transcript exons" ] ++ map (BS.unpack . repr) clocs

cdsLoc :: SpliceSeqLoc -> [ContigSeqLoc] -> Either String (Maybe Loc.ContigLoc)                
cdsLoc _ [] = return Nothing
cdsLoc (OnSeq trxname trxloc) cdses@(_:_) = do
  (seqname, rawcontigs) <- allSameName cdses
  when (trxname /= seqname) $ Left . unwords $ [ "CDS sequence name mismatch", 
                                                 show $ unSeqLabel trxname
                                               , show $ unSeqLabel seqname ]
  contigs <- sortclocs rawcontigs
  (contigIn0:contigInRest) <- mapM (cdsIntoTranscript trxloc) contigs
  cloc <- foldM mergeCLocs contigIn0 contigInRest
  return $! Just $! Loc.extend (0, 3) cloc -- Include the stop codon
  where sortclocs clocs = maybe badClocs Right . sortContigs $ clocs
          where badClocs = Left . unwords $ [ "bad transcript CDSes" ] ++ map (BS.unpack . repr) clocs  

cdsIntoTranscript :: SpLoc.SpliceLoc -> Loc.ContigLoc -> Either String Loc.ContigLoc
cdsIntoTranscript trxloc cdscontig = maybe badIn Right . Loc.clocInto cdscontig $ trxloc
  where badIn = Left . unwords $ [ "Mapping CDS contig into transcript: "
                                 , BS.unpack . repr $ cdscontig
                                 , BS.unpack . repr $ trxloc
                                 ]

mergeCLocs :: Loc.ContigLoc -> Loc.ContigLoc -> Either String Loc.ContigLoc
mergeCLocs cloc0 clocnext
  | (Loc.strand cloc0 == Plus) && (Loc.startPos clocnext == Loc.endPos cloc0 `Pos.slide` 1)
    = return $! Loc.extend (0, Loc.length clocnext) cloc0
  | otherwise = Left . unwords $ [ "Merging non-adjacent contigs: "
                                 , BS.unpack . repr $ cloc0
                                 , BS.unpack . repr $ clocnext
                                 ]
                                   

allSameName :: [OnSeq a] -> Either String (SeqLabel, [a])
allSameName s = case group . map onSeqLabel $ s of
  [(name0:_)] -> return (name0, map unOnSeq s)
  names -> Left $ "allSameName: names " ++ show (map (unSeqLabel . head) names)
  
data GtfLine = GtfLine { gtfGene, gtfTrx, gtfFtype :: !BS.ByteString, gtfLoc :: !ContigSeqLoc }

data GtfTrxs = GtfTrxs { gtfExonLocs, gtfCdsLocs :: !(HM.HashMap BS.ByteString [ContigSeqLoc])
                       , gtfTogene :: !(HM.HashMap BS.ByteString BS.ByteString)
                       }
data GtfAttr = GtfAttr { attrName, attrValue :: !BS.ByteString } deriving (Show)
               
ftypeCds :: BS.ByteString
ftypeCds = BS.pack "CDS"

ftypeExon :: BS.ByteString
ftypeExon = BS.pack "exon"

insertGtfLine :: GtfTrxs -> GtfLine -> GtfTrxs
insertGtfLine trxsin l 
  | gtfFtype l == ftypeCds  = insertCds  . insertGene $ trxsin
  | gtfFtype l == ftypeExon = insertExon . insertGene $ trxsin
  | otherwise = trxsin
    where trx = BS.copy $ gtfTrx l
          insertGene t0 = {-# SCC "insertGene" #-} t0 { gtfTogene = HM.insert trx (BS.copy . gtfGene $ l) (gtfTogene t0) }
          insertCds t0 = {-# SCC "insertCds" #-} t0 { gtfCdsLocs = HM.insertWith (++) trx [gtfLoc l] (gtfCdsLocs t0) }
          insertExon t0 = {-# SCC "insertExon" #-} t0 { gtfExonLocs = HM.insertWith (++) trx [gtfLoc l] (gtfExonLocs t0) } 

gtflineIter :: (Monad m) => Iter.Enumeratee [BS.ByteString] [GtfLine] m a
gtflineIter = Iter.convStream $ Iter.head >>= liftM (: []) . handleErr . gtfline
  where handleErr = either (Iter.throwErr . Iter.iterStrExc) return 

gtfline :: BS.ByteString -> Either String GtfLine
gtfline l = addErrorLine $ case V.fromList . BS.split '\t' $ l of
  fields | V.length fields < 9 -> Left $ "Expected 9 fields, saw " ++ show (V.length fields)
         | otherwise -> let !seqname = toSeqLabel . BS.copy $ fields V.! 0
                            !ftype = BS.copy $ fields V.! 2
                        in do start <- decode $ fields V.! 3
                              end <- decode $ fields V.! 4
                              str <- case fields V.! 6 of { "+" -> Right Plus; "-" -> Right Minus; ch -> Left ("Bad strand character " ++ show ch) }
                              let !loc = Loc.fromBoundsStrand (start - 1) (end - 1) str
                                  !seqloc = OnSeq seqname loc
                              attrs <- gtfattrs $ fields V.! 8
                              gene <- reqattr "gene_id" attrs
                              trx <- reqattr "transcript_id" attrs
                              return $! GtfLine gene trx ftype seqloc
  where addErrorLine good@(Right _) = good
        addErrorLine (Left err) = Left $! err ++ "\nParsing GTF line " ++ show l
        reqattr name attrs = case find ((== name) . attrName) attrs of
          Just (GtfAttr _name value) -> Right value
          Nothing -> Left $ "Missing required attribute " ++ show name ++ " within " ++ show (map attrName attrs)

gtfattrs :: BS.ByteString -> Either String [GtfAttr]
gtfattrs attrstr = C.unfoldM gtfattr attrstr C.$$ C.consume

gtfattr :: BS.ByteString -> Either String (Maybe (GtfAttr, BS.ByteString))
gtfattr l0 = case BS.dropWhile isSpace l0 of
  l | BS.null l -> return Nothing
    | otherwise -> let !(name0, valrest) = BS.break isSpace l
                       !name = BS.copy name0
                   in do (value0, rest) <- gtfvalue $ BS.dropWhile isSpace valrest
                         let !value = BS.copy value0
                             !attr = GtfAttr name value
                         return $! Just (attr, rest)
  where gtfvalue valrest = case BS.uncons valrest of
          Just ('\"', str) -> let (val, rest0) = BS.break (== '\"') str
                              in if BS.isPrefixOf "\";" rest0
                                 then return (val, BS.drop 2 rest0)
                                 else Left $ "Malformed quoted attribute value starting at\n" ++ show valrest
          _ -> let (val, rest0) = BS.break (== ';') valrest
               in if BS.isPrefixOf ";" rest0
                  then return (val, BS.drop 1 rest0)
                  else Left $ "Malformed unquoted attribute value starting at\n" ++ show valrest
