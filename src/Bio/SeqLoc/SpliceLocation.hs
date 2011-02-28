{-# LANGUAGE TypeFamilies, FlexibleContexts #-}

{-| Data type for a more general sequence location consiting of
disjoint ranges of positions on a sequence.

Throughout, /sequence position/ refers to a 'Pos.Pos' which includes a
strand, as opposed to an /offset/, which refers to a 'Pos.Offset' with
no strand.
 -}

module Bio.SeqLoc.SpliceLocation ( 
  -- * Sequence locations
  SpliceLoc(..)
  , fromContigs
  ) where 

import Prelude hiding (length)

import Control.Applicative
import Control.Arrow ((***))
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List (foldl')
import Data.Maybe

import qualified Data.Attoparsec.Char8 as AP

import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.Location
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import qualified Bio.SeqLoc.SeqData as SeqData

-- | General (disjoint) sequence region consisting of a concatenated
-- set of one or more contiguous regions.
data SpliceLoc = SpliceLocLast { contig :: !ContigLoc } 
               | SpliceLocPrev { contig :: !ContigLoc
                               , next :: !SpliceLoc
                               }
               deriving (Eq, Ord, Show)

fromContigs :: [ContigLoc] -> SpliceLoc
fromContigs [] = error $ "Bio.SeqLoc.SpliceLocation.fromContigs: empty contigs"
fromContigs [l] = SpliceLocLast l
fromContigs (l:rest) = SpliceLocPrev l (fromContigs rest)

tails :: SpliceLoc -> [SpliceLoc]
tails sll@(SpliceLocLast _) = [sll]
tails slp@(SpliceLocPrev _ n) = slp : (tails n)

instance Stranded SpliceLoc where
  revCompl sll = foldl' addprev newlast slrest
    where (sl0:slrest) = tails sll
          newlast = SpliceLocLast (revCompl . contig $ sl0)
          addprev slnew slold = SpliceLocPrev (revCompl . contig $ slold) slnew

instance LocRepr SpliceLoc where
  repr = BS.intercalate (BS.singleton ';') . map repr . contigs
  unrepr = fromContigs <$> AP.sepBy1 unrepr (AP.char ';')

instance Location SpliceLoc where
  length = foldl' (\len c -> len + length c) 0 . contigs
  bounds = (minimum *** maximum) . unzip . map bounds . contigs  
  seqData sequ = liftM SeqData.concat . mapM (seqData sequ) . contigs
  seqDataPad sequ = SeqData.concat . map (seqDataPad sequ) . contigs
  posInto pos = posIntoContigs pos . contigs
  posOutof pos = posOutofContigs pos . contigs
  startPos = startPos . firstContig
  endPos = endPos . lastContig
  clocInto = slocClocInto
  clocOutof = slocClocOutof
  extend = slocExtend
  offsetWithin off = or . map (offsetWithin off) . contigs
  posWithin pos = or . map (posWithin pos) . contigs
  contigOverlaps c = any (contigOverlaps c ) . contigs
  toContigs = contigs

contigs :: SpliceLoc -> [ContigLoc]
contigs (SpliceLocLast c) = [c]
contigs (SpliceLocPrev c n) = c : contigs n

firstContig :: SpliceLoc -> ContigLoc
firstContig = contig

lastContig :: SpliceLoc -> ContigLoc
lastContig (SpliceLocLast c) = c
lastContig (SpliceLocPrev _ n) = lastContig n

--reloffsets :: SpliceLoc -> [Pos.Offset]
--reloffsets = init . scanl (\off c -> off + length c) 0 . contigs

contigsAndOffsets :: SpliceLoc -> [(ContigLoc, Pos.Offset)]
contigsAndOffsets = go 0
  where go off (SpliceLocLast c) = [(c, off)]
        go off (SpliceLocPrev c n) = (c, off) : (go (off + length c) n)

posIntoContigs :: Pos.Pos -> [ContigLoc] -> Maybe Pos.Pos
posIntoContigs pos = go 0
  where go _ [] = Nothing
        go dlen (c:rest) = maybe onRest onContig . posInto pos $ c
          where onContig pin = Just $! pin `Pos.slide` dlen
                onRest = go (dlen + length c) rest

posOutofContigs :: Pos.Pos -> [ContigLoc] -> Maybe Pos.Pos
posOutofContigs _ [] = Nothing
posOutofContigs pos (c:rest) = maybe onRest onContig . posOutof pos $ c
  where onRest = posOutofContigs (Pos.slide pos . negate . length $ c) rest
        onContig = Just

slocClocInto :: ContigLoc -> SpliceLoc -> Maybe ContigLoc
slocClocInto subcloc = listToMaybe . catMaybes . map into . contigsAndOffsets
    where into (cloc, off) = liftM (slide off) . clocInto subcloc $ cloc

slocClocOutof :: ContigLoc -> SpliceLoc -> Maybe SpliceLoc
slocClocOutof (ContigLoc pos5 len cstrand) 
    = liftM (stranded cstrand . fromContigs) . regionOutof pos5 len . contigs

regionOutof :: Pos.Offset -> Pos.Offset -> [ContigLoc] -> Maybe [ContigLoc]
regionOutof _ _ [] = Nothing
regionOutof pos5 len (cloc0:rest)
    | pos5 < 0 = Nothing
    | len < 0 = error $ "Bio.SeqLoc.SpliceLocation.regionOutof: input, " ++ show (pos5, len, cloc0)
    | pos5 >= len0 = regionOutof (pos5 - len0) len rest
    | (pos5 + len <= len0) = case clocOutof (ContigLoc pos5 len Fwd) cloc0 of
                               Just out0 -> Just [out0]
                               Nothing -> error $ "regionOutof: final bounds failure, " ++ show (pos5, len, cloc0)
    | otherwise = let outlen0 = len0 - pos5
                  in case clocOutof (ContigLoc pos5 outlen0 Fwd) cloc0 of
                       Just out0 -> liftM (out0 :) . regionOutof 0 (len - outlen0) $ rest
                       Nothing -> error $ "regionOutof: internal bounds failure, " ++ show (pos5, len, outlen0, cloc0)
    where len0 = length cloc0

slocExtend :: (Pos.Offset, Pos.Offset) -> SpliceLoc -> SpliceLoc
slocExtend (ext5, ext3) sloc = extended3 { contig = extend (ext5, 0) . contig $ extended3 }
  where extend3 (SpliceLocPrev c n) = SpliceLocPrev c (extend3 n)  
        extend3 (SpliceLocLast c) = SpliceLocLast $ extend (0, ext3) c
        extended3 = extend3 sloc
