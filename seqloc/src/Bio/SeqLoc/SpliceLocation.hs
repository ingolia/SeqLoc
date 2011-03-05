{-# LANGUAGE TypeFamilies, FlexibleContexts, OverloadedStrings #-}

{-| Data type for a more general sequence location consiting of
disjoint ranges of positions on a sequence.

Throughout, /sequence position/ refers to a 'Pos.Pos' which includes a
strand, as opposed to an /offset/, which refers to a 'Pos.Offset' with
no strand.
 -}

module Bio.SeqLoc.SpliceLocation ( 
  -- * Sequence locations
  SpliceLoc
  , fromContigs
  , locOutof, locInto
  , mergeContigs, mergeAdjContigs
  ) where 

import Prelude hiding (length)

import Control.Applicative
import Control.Arrow ((***))
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List (foldl')
import Data.Maybe

import qualified Data.Attoparsec.Zepto as ZP

import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.Location
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import qualified Bio.SeqLoc.SeqData as SeqData

-- | General (disjoint) sequence region consisting of a concatenated
-- set of one or more contiguous regions.
data SpliceLoc = SpliceLocLast { contig :: !ContigLoc } 
               | SpliceLocPrev { contig :: !ContigLoc
                               , _next :: !SpliceLoc
                               }
               deriving (Eq, Ord, Show)

singleton :: ContigLoc -> SpliceLoc
singleton = SpliceLocLast

consContig :: ContigLoc -> SpliceLoc -> Maybe SpliceLoc
consContig c sploc | goodJunction c (contig sploc) = Just $! SpliceLocPrev c sploc
                   | otherwise = Nothing
                     
goodJunction :: ContigLoc -> ContigLoc -> Bool
goodJunction c5 c3 = sameStrand && inOrder
  where sameStrand = strand c5 == strand c3
        c5end = snd . bounds $ c5
        c3start = fst . bounds $ c3
        inOrder = case strand c5 of
          Fwd -> c5end < c3start
          RevCompl -> c5end > c3start

fromContigs :: [ContigLoc] -> Maybe SpliceLoc
fromContigs [] = Nothing
fromContigs [l] = Just $! singleton l
fromContigs (l:rest) = fromContigs rest >>= consContig l

tails :: SpliceLoc -> [SpliceLoc]
tails sll@(SpliceLocLast _) = [sll]
tails slp@(SpliceLocPrev _ n) = slp : (tails n)

instance Stranded SpliceLoc where
  revCompl sll = fromMaybe badRevCompl . foldl' addprev newlast $ slrest
    where (sl0:slrest) = tails sll
          newlast = Just $! singleton . revCompl . contig $ sl0
          addprev slnew slold = slnew >>= consContig (revCompl . contig $ slold)
          badRevCompl = error $ "Bad junction doing reverse complement on " ++ (BS.unpack . repr) sll

instance LocRepr SpliceLoc where
  repr = BS.intercalate (BS.singleton ';') . map repr . contigs
  unrepr = (fromContigs <$> scan) >>= maybe (fail "bad contig order") return
    where scan = liftA2 (:) unrepr ((ZP.string ";" *> scan) <|> pure [])

instance Location SpliceLoc where
  strand = strand . contig
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

locOutof :: (Location l) => SpliceLoc -> l -> Maybe SpliceLoc
locOutof sploc outer = mapM (flip clocOutof outer) (toContigs sploc) >>= 
                       fromContigs . concat . map toContigs

locInto :: (Location l) => SpliceLoc -> l -> Maybe SpliceLoc
locInto sploc outer = mapM (flip clocInto outer) (toContigs sploc) >>=
                      fromContigs . mergeContigs . concat . map toContigs

mergeContigs :: [ContigLoc] -> [ContigLoc]
mergeContigs [] = []
mergeContigs [clast] = [clast]
mergeContigs (cprev:rest@(_:_)) = case mergeContigs rest of
  [] -> error $ "mergeContigs: empty rest"
  mergerest@(cnext:afternext) -> case mergeAdjContigs cprev cnext of
    Just cmerge -> cmerge : afternext
    Nothing     -> cprev : mergerest

mergeAdjContigs :: ContigLoc -> ContigLoc -> Maybe ContigLoc
mergeAdjContigs clocprev clocnext 
  | startPos clocnext == endPos (extend (0, 1) clocprev)
     = Just $! fromPosLen (startPos clocprev) (length clocprev + length clocnext)
  | otherwise = Nothing

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
slocClocOutof cloc
    = liftM (stranded (strand cloc) . fromContigsErr) . regionOutof (offset5 cloc) (length cloc) . contigs
  where fromContigsErr ctgs = fromMaybe badContigs $! fromContigs ctgs
          where badContigs = error . unwords $ 
                             [ "bad contig order in slocClocOutof" ] ++
                             map (BS.unpack . repr) ctgs
        
  
regionOutof :: Pos.Offset -> Pos.Offset -> [ContigLoc] -> Maybe [ContigLoc]
regionOutof _ _ [] = Nothing
regionOutof pos5 len (cloc0:rest)
    | pos5 < 0 = Nothing
    | len < 0 = error $ "Bio.SeqLoc.SpliceLocation.regionOutof: input, " ++ show (pos5, len, cloc0)
    | pos5 >= len0 = regionOutof (pos5 - len0) len rest
    | (pos5 + len <= len0) = let subcloc = fromPosLen (Pos.Pos pos5 Fwd) len
                             in case clocOutof subcloc cloc0 of
                               Just out0 -> Just [out0]
                               Nothing -> error $ "regionOutof: final bounds failure, " ++ (BS.unpack . repr) subcloc
    | otherwise = let outlen0 = len0 - pos5
                      subcloc = fromPosLen (Pos.Pos pos5 Fwd) outlen0
                  in case clocOutof subcloc cloc0 of
                       Just out0 -> liftM (out0 :) . regionOutof 0 (len - outlen0) $ rest
                       Nothing -> error $ "regionOutof: internal bounds failure, " ++ (BS.unpack . repr) subcloc
    where len0 = length cloc0

slocExtend :: (Pos.Offset, Pos.Offset) -> SpliceLoc -> SpliceLoc
slocExtend (ext5, ext3) sloc = extended3 { contig = extend (ext5, 0) . contig $ extended3 }
  where extend3 (SpliceLocPrev c n) = SpliceLocPrev c (extend3 n)  
        extend3 (SpliceLocLast c) = SpliceLocLast $ extend (0, ext3) c
        extended3 = extend3 sloc
