{-# LANGUAGE BangPatterns, MagicHash, OverloadedStrings, UnboxedTuples #-}

module Bio.SeqLoc.GTF
       ( readGtfTranscripts
       )
       where 

import Control.Applicative
import Control.Monad
import qualified Data.ByteString as BSW
import qualified Data.ByteString.Char8 as BS
import Data.ByteString.Internal (c2w)
import Data.List
import Data.Maybe

import qualified Data.Attoparsec.Char8 as AP (isDigit_w8, isSpace_w8)
--import qualified Data.Attoparsec.Char8 as AP
import qualified Data.Attoparsec.Zepto as ZP
import qualified Data.HashMap.Strict as HM
import qualified Data.Iteratee as Iter
import qualified Data.Iteratee.Char as IterChar
--import qualified Data.Iteratee.Exception as Iter

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

readGtfTranscripts :: FilePath -> IO [Transcript]
readGtfTranscripts = Iter.fileDriver gtfTrxsIter >=> 
                     either (ioError . userError) return . mkTranscripts
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

mkTranscript :: BS.ByteString -> [ContigSeqLoc] -> [ContigSeqLoc] -> BS.ByteString -> Either String Transcript
mkTranscript trx exons cdses gene = moderr $ do loc <- exonLoc exons
                                                cdsloc <- cdsLoc loc cdses
                                                return $! Transcript (SeqName gene) (SeqName trx) loc cdsloc
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
  when (trxname /= seqname) $ Left . unwords $ [ "CDS sequence name mismatch", show trxname, show seqname ]
  contigs <- sortclocs rawcontigs
  (contigIn0:contigInRest) <- mapM (cdsIntoTranscript trxloc) contigs
  liftM Just $! foldM mergeCLocs contigIn0 contigInRest
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
  | (Loc.strand cloc0 == Fwd) && (Loc.startPos clocnext == Loc.endPos cloc0 `Pos.slide` 1)
    = return $! Loc.extend (0, Loc.length clocnext) cloc0
  | otherwise = Left . unwords $ [ "Merging non-adjacent contigs: "
                                 , BS.unpack . repr $ cloc0
                                 , BS.unpack . repr $ clocnext
                                 ]
                                   

allSameName :: [OnSeq a] -> Either String (SeqName, [a])
allSameName s = case group . map onSeqName $ s of
  [(name0:_)] -> return (name0, map unOnSeq s)
  names -> Left $ "allSameName: names " ++ show (map (unSeqName . head) names)
  
data GtfLine = GtfLine { gtfGene, gtfTrx, gtfFtype :: !BS.ByteString, gtfLoc :: !ContigSeqLoc } deriving (Show)

data GtfTrxs = GtfTrxs { gtfExonLocs, gtfCdsLocs :: !(HM.HashMap BS.ByteString [ContigSeqLoc])
                       , gtfTogene :: !(HM.HashMap BS.ByteString BS.ByteString)
                       } deriving (Show)
               
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
gtflineIter = Iter.convStream $ Iter.head >>= liftM (: []) . handleErr . ZP.parse gtfline
  where handleErr = either (Iter.throwErr . Iter.iterStrExc) return 

-- Does NOT consume the remainder of the line
gtfline :: ZP.Parser GtfLine
gtfline = do seqname <- field
             _source <- dropField
             ftype <- field
             start <- decfield
             end <- decfield
             _score <- dropField
             str <- strand
             _frame <- dropField
             gene <- attr "gene_id" <* ZP.string (BS.pack "; ")
             trx <- attr "transcript_id" <* ZP.string (BS.pack "; ")
             let name = SeqName . BS.copy $ seqname
                 loc = Loc.fromBoundsStrand (start - 1) (end - 1) str
             return $! GtfLine gene trx ftype (OnSeq name loc)

strand :: ZP.Parser Strand
strand = ((ZP.string "+\t" *> return Fwd) <|>
          (ZP.string "-\t" *> return RevCompl))

attr :: String -> ZP.Parser BS.ByteString
attr name = ZP.string (BS.pack name) *> ZP.takeWhile AP.isSpace_w8 *> ZP.string "\"" *>
            ZP.takeWhile (/= c2w '\"') <* ZP.string "\""

decfield :: (Integral a) => ZP.Parser a
decfield = decimal <* ZP.string "\t"

field :: ZP.Parser BS.ByteString
field = ZP.takeWhile (/= c2w '\t') <* ZP.string "\t"

dropField :: ZP.Parser ()
dropField = ZP.takeWhile (/= c2w '\t') *> ZP.string "\t" *> return ()

decimal :: (Integral a) => ZP.Parser a
decimal = decode <$> ZP.takeWhile (AP.isDigit_w8)
  where decode = fromIntegral . BSW.foldl' step (0 :: Int)
        step a w = a * 10 + fromIntegral (w - 48)
                           
