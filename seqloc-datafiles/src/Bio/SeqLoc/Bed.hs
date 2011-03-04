{-# LANGUAGE OverloadedStrings #-}

module Bio.SeqLoc.Bed
       ( readBedTranscripts
       , bedZP, bedTranscriptEnum
       , transcriptToBed, transcriptToBedStd
       )
       where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Ord

import qualified Data.Attoparsec.Zepto as ZP
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

transcriptToBedStd :: Transcript -> BS.ByteString
transcriptToBedStd = transcriptToBed "0" "0"

transcriptToBed :: BS.ByteString -> BS.ByteString -> Transcript -> BS.ByteString
transcriptToBed score rgb trx = unfields fields
  where unfields = BS.intercalate (BS.singleton '\t')
        fields = [ unSeqName chrom
                 , repr $ chromStart
                 , repr $ chromEnd + 1
                 , unSeqName . trxId $ trx
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
        strandchr = case Loc.strand loc of Fwd -> "+"; RevCompl -> "-"
        (thickStart, thickEnd) = maybe noCds (Loc.bounds . unOnSeq) . cdsLocation $ trx
        noCds = (chromStart, chromStart - 1)
        contigs = sortBy (comparing Loc.offset5) . Loc.toContigs $ loc
        blockSizes = map Loc.length contigs
        blockStarts = map (subtract chromStart . Loc.offset5) contigs
        unCommaList = BS.concat . map (flip BS.append (BS.singleton ',') . repr)

readBedTranscripts :: FilePath -> IO [Transcript]
readBedTranscripts = Iter.fileDriver (bedTranscriptEnum Iter.stream2list)
                     
bedTranscriptEnum :: (Monad m) => Iter.Iteratee [Transcript] m a -> Iter.Iteratee BS.ByteString m a
bedTranscriptEnum = Iter.joinI . IterChar.enumLinesBS . Iter.joinI . bedLineEnum

bedLineEnum :: (Monad m) => Iter.Enumeratee [BS.ByteString] [Transcript] m a
bedLineEnum = Iter.convStream $ Iter.head >>= liftM (: []) . handleErr . ZP.parse bedZP
  where handleErr = either (Iter.throwErr . Iter.iterStrExc) return 

bedZP :: ZP.Parser Transcript
bedZP = do chrom <- field -- The name of the chromosome
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
           _score <- dropField -- A score between 0 and 1000.
           str <- strand -- Defines the strand
           thickStart <- decfield -- The starting position at which
                                  -- the feature is drawn thickly (for
                                  -- example, the start codon in gene
                                  -- displays).
           thickEnd <- decfield -- The ending position at which the
                                -- feature is drawn thickly (for
                                -- example, the stop codon in gene
                                -- displays).
           _itemRGB <- dropField -- An RGB value of the form R,G,B
                                 -- (e.g. 255,0,0).
           blockCount <- decfield -- The number of blocks (exons) in
                                  -- the BED line.
           blockSizes <- commas blockCount decimal <* dropField
                         -- A comma-separated list of the block sizes.                         
           blockStarts <- commas blockCount decimal 
                          -- A comma-separated list of block starts.
           loc <- bedTrxLoc chromStart chromEnd str $ zip blockSizes blockStarts
           unless (Loc.bounds loc == (chromStart, chromEnd - 1)) $
             fail $ "Bio.SeqLoc.Bed: bad sploc:" ++ 
             (BS.unpack . BS.unwords $ [ repr loc, repr chromStart, repr chromEnd ])
           cdsloc <- if thickStart >= thickEnd
                        then return Nothing
                        else liftM Just $! bedCdsLoc loc thickStart thickEnd
           let n = SeqName $ BS.copy name
               c = SeqName $ BS.copy chrom
           return $! Transcript n n (OnSeq c loc) cdsloc
           
bedTrxLoc :: (Monad m) => Pos.Offset -> Pos.Offset -> Strand -> [(Pos.Offset, Pos.Offset)] -> m SpLoc.SpliceLoc
bedTrxLoc chromStart chromEnd str = maybe badContigs (return . stranded str) . 
                                    SpLoc.fromContigs . map blockContig
  where blockContig (bsize, bstart) = Loc.fromPosLen (Pos.Pos (chromStart + bstart) Fwd) bsize
        badContigs = fail $ "Bio.SeqLoc.Bed: bad blocks in " ++ show (chromStart, chromEnd)
        
bedCdsLoc :: (Monad m) => SpLoc.SpliceLoc -> Pos.Offset -> Pos.Offset -> m Loc.ContigLoc
bedCdsLoc loc thickStart thickEnd 
  = maybe badCdsLoc return $ do
    relstart <- Loc.posInto (Pos.Pos thickStart Fwd) loc
    relend <- Loc.posInto (Pos.Pos (thickEnd - 1) Fwd) loc
    return $! stranded (Loc.strand loc) $ Loc.fromStartEnd (Pos.offset relstart) (Pos.offset relend)
      where badCdsLoc = fail $ "Bio.SeqLoc.Bed: bad cds in " ++ 
                        (BS.unpack . BS.unwords $ [ repr loc, repr thickStart, repr thickEnd ])

commas :: Int -> ZP.Parser a -> ZP.Parser [a]
commas n p | n < 1     = return [] 
           | otherwise = (:) <$> p <*>
                         replicateM (n - 1) (ZP.string "," *> p)