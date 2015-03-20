{-# LANGUAGE OverloadedStrings, FlexibleContexts #-}

{-| Utilities for reading and writing BED format gene annotations -}
module Bio.SeqLoc.Bed
       ( readBedTranscripts
       , bedZP, bedTranscriptEnum
       , bedConduit, unbedConduit
       , transcriptToBed, transcriptToBedStd
       )
       where

import Control.Applicative
import qualified Control.Exception.Lifted as E
import Control.Monad
import Control.Monad.Base
import Control.Monad.Trans.Resource
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import Data.Ord

import qualified Data.Attoparsec.Zepto as ZP
import qualified Data.Conduit as C
import qualified Data.Conduit.Binary as CB
import qualified Data.Conduit.List as C
import qualified Data.Iteratee as Iter
import qualified Data.Iteratee.Char as IterChar

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import Bio.SeqLoc.ZeptoUtils

-- | Convert a 'Transcript' to a BED annotation line.
transcriptToBedStd :: Transcript -> BS.ByteString
transcriptToBedStd = transcriptToBed "0" "0"

-- | Convert a 'Transcript' to a BED annotation line, specifying the
-- /score/ and /itemRGB/ fields.
transcriptToBed :: BS.ByteString -- ^ score
                   -> BS.ByteString -- ^ itemRGB
                   -> Transcript -- ^ transcript
                   -> BS.ByteString
transcriptToBed score rgb trx = unfields fields
  where unfields = BS.intercalate (BS.singleton '\t')
        fields = [ unSeqLabel chrom
                 , repr $ chromStart
                 , repr $ chromEnd + 1
                 , unSeqLabel . trxId $ trx
                 , score
                 , strandchr
                 , repr $ thickStart
                 , repr $ thickEnd + 1
                 , rgb
                 , BS.pack . show . length $ blockSizes
                 , unCommaList blockSizes
                 , unCommaList blockStarts
                 ]
        (OnSeq chrom loc) = location trx
        (chromStart, chromEnd) = Loc.bounds loc
        strandchr = case Loc.strand loc of Plus -> "+"; Minus -> "-"
        (thickStart, thickEnd) = maybe noCds (Loc.bounds . unOnSeq) . cdsLocation $ trx
        noCds = (chromStart, chromStart - 1)
        contigs = sortBy (comparing Loc.offset5) . Loc.toContigs $ loc
        blockSizes = map Loc.length contigs
        blockStarts = map (subtract chromStart . Loc.offset5) contigs
        unCommaList = BS.concat . map (flip BS.append (BS.singleton ',') . repr)

-- | Read all BED format annotations in a BED file
readBedTranscripts :: FilePath -> IO [Transcript]
readBedTranscripts bedfile = runResourceT $ C.runConduit $
                             CB.sourceFile bedfile C.$= bedConduit C.$= C.consume
                      
-- | Iteratee to convert an 'Iter.Iteratee' over a 'BS.ByteString',
-- such as the standard 'Iter.fileDriver', into an iteratee over a
-- list of 'Transcript' annotations from the file.
bedTranscriptEnum :: (Monad m) => Iter.Iteratee [Transcript] m a -> Iter.Iteratee BS.ByteString m a
bedTranscriptEnum = Iter.joinI . IterChar.enumLinesBS . Iter.joinI . bedLineEnum

bedLineEnum :: (Monad m) => Iter.Enumeratee [BS.ByteString] [Transcript] m a
bedLineEnum = Iter.convStream $ Iter.head >>= liftM (: []) . handleErr . ZP.parse bedZP
  where handleErr = either (Iter.throwErr . Iter.iterStrExc) return 

-- | Conduit from a 'BS.ByteString' source such as a BED file to a
-- source of 'Transcript' annotations from the file.
bedConduit :: (Monad m, MonadBase IO m) => C.Conduit BS.ByteString m Transcript
bedConduit = CB.lines C.$= loop
  where loop = C.head >>= maybe (return ())
               (\l -> case ZP.parse bedZP l of
                   Left err -> E.ioError . userError $ err ++ "\n  in BED line\n"  ++ show l
                   Right res -> C.yield res >> loop)

unbedConduit :: (Monad m) => C.Conduit Transcript m BS.ByteString
unbedConduit = C.head >>= maybe (return ()) go
  where go t = C.yield (transcriptToBedStd t `BS.append` "\n") >> unbedConduit

-- | Minimalistic 'ZP.Parser'-style parser for a BED format line, not
-- including the trailing newline.
bedZP :: ZP.Parser Transcript
bedZP = do chrom <- firstfield -- The name of the chromosome
           chromStart <- decfield -- The starting position of the
                                  -- feature in the chromosome or
                                  -- scaffold. The first base in a
                                  -- chromosome is numbered 0
           chromEnd <- decfield -- The ending position of the feature
                                -- in the chromosome or scaffold. The
                                -- chromEnd base is not included in
                                -- the display of the feature. For
                                -- example, the first 100 bases of a
                                -- chromosome are defined as
                                -- chromStart=0, chromEnd=100, and
                                -- span the bases numbered 0-99.
           name <- field -- Defines the name of the BED line.
           _score <- unlessAtEnd dropField -- A score between 0 and 1000.
           str <- fromMaybe Plus <$> unlessAtEnd strand -- Defines the strand
           thickStart <- unlessAtEnd decfield -- The starting position at which
                                  -- the feature is drawn thickly (for
                                  -- example, the start codon in gene
                                  -- displays).
           thickEnd <- unlessAtEnd decfield -- The ending position at which the
                                -- feature is drawn thickly (for
                                -- example, the stop codon in gene
                                -- displays).
           _itemRGB <- unlessAtEnd dropField -- An RGB value of the form R,G,B
                                             -- (e.g. 255,0,0).
           blockCount <- unlessAtEnd decfield -- The number of blocks (exons) in
                                  -- the BED line.
           blockSizes <- case blockCount of
             Just ct -> unlessAtEnd $ commas ct decimal
                        -- A comma-separated list of the block sizes.                         
             Nothing -> return Nothing
           blockStarts <- case blockCount of
             Just ct -> unlessAtEnd $ commas ct decimal 
                        -- A comma-separated list of block starts.
             Nothing -> return Nothing
           loc <- case liftM2 zip blockSizes blockStarts of
             Just blocks -> bedTrxLoc chromStart chromEnd str blocks
             Nothing -> maybe badContig return $ 
                        SpLoc.fromContigs [ Loc.fromBoundsStrand chromStart (chromEnd - 1) str ]
                          where badContig = error "bedZP: Bad singleton sploc!"
           unless (Loc.bounds loc == (chromStart, chromEnd - 1)) $
             fail $ "Bio.SeqLoc.Bed: bad sploc:" ++ 
             (BS.unpack . BS.unwords $ [ repr loc, repr chromStart, repr chromEnd ])
           cdsloc <- case liftM2 (,) thickStart thickEnd of
             Just (start, end) | end > start -> liftM Just $! bedCdsLoc loc start end
             _ -> return Nothing
           let n = toSeqLabel $ BS.copy name
               c = toSeqLabel $ BS.copy chrom
           return $! Transcript n n (OnSeq c loc) cdsloc
           
bedTrxLoc :: (Monad m) => Pos.Offset -> Pos.Offset -> Strand -> [(Pos.Offset, Pos.Offset)] -> m SpLoc.SpliceLoc
bedTrxLoc chromStart chromEnd str = maybe badContigs (return . stranded str) . 
                                    SpLoc.fromContigs . map blockContig
  where blockContig (bsize, bstart) = Loc.fromPosLen (Pos.Pos (chromStart + bstart) Plus) bsize
        badContigs = fail $ "Bio.SeqLoc.Bed: bad blocks in " ++ show (chromStart, chromEnd)
        
bedCdsLoc :: (Monad m) => SpLoc.SpliceLoc -> Pos.Offset -> Pos.Offset -> m Loc.ContigLoc
bedCdsLoc loc thickStart thickEnd 
  = maybe badCdsLoc return $ do
    relstart <- Loc.posInto (Pos.Pos thickStart Plus) loc
    relend <- Loc.posInto (Pos.Pos (thickEnd - 1) Plus) loc
    case Loc.strand loc of
      Plus -> if relstart <= relend
                 then return $! Loc.fromBoundsStrand (Pos.offset relstart) (Pos.offset relend) Plus
                 else Nothing
      Minus -> if relend <= relstart
                  then return $! Loc.fromBoundsStrand (Pos.offset relend) (Pos.offset relstart) Plus
                  else Nothing
      where badCdsLoc = fail $ "Bio.SeqLoc.Bed: bad cds in " ++ 
                        (BS.unpack . BS.unwords $ [ repr loc, repr thickStart, repr thickEnd ])

commas :: Int -> ZP.Parser a -> ZP.Parser [a]
commas n p | n < 1     = return [] 
           | otherwise = ZP.string "\t" *> 
                         ( (:) <$> p <*>
                           replicateM (n - 1) (ZP.string "," *> p) )
                         <* (ZP.string "," <|> return ())
