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

  -- * Locations and positions
  , bounds, length, startPos, endPos -- , posInto, posOutof
--  , clocInto, clocOutof, isWithin, overlaps

  -- * Extracting subsequences
--  , seqDataM, seqDataPadded

  -- * Transforming locations
--  , extend

  -- * Displaying locations
--  , parser
  ) where 

import Prelude hiding (length)

import Control.Applicative
import Control.Arrow ((***))
import Data.List (foldl', mapAccumL)
import Data.Maybe

import qualified Bio.SeqLoc.ContigLocation as CLoc
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import qualified Bio.SeqLoc.SeqData as SeqData

-- | General (disjoint) sequence region consisting of a concatenated
-- set of one or more contiguous regions.
data SpliceLoc = SpliceLocLast { contig :: !CLoc.ContigLoc } 
               | SpliceLocPrev { contig :: !CLoc.ContigLoc
                               , next :: !SpliceLoc
                               }
               deriving (Eq, Ord, Show)

tails :: SpliceLoc -> [SpliceLoc]
tails sll@(SpliceLocLast _) = [sll]
tails slp@(SpliceLocPrev _ n) = slp : (tails n)

instance Stranded SpliceLoc where
  revCompl sll = foldl' addprev newlast slrest
    where (sl0:slrest) = tails sll
          newlast = SpliceLocLast (revCompl . contig $ sl0)
          addprev slnew slold = SpliceLocPrev (revCompl . contig $ slold) slnew

contigs :: SpliceLoc -> [CLoc.ContigLoc]
contigs (SpliceLocLast c) = [c]
contigs (SpliceLocPrev c n) = c : contigs n

firstContig :: SpliceLoc -> CLoc.ContigLoc
firstContig = contig

lastContig :: SpliceLoc -> CLoc.ContigLoc
lastContig (SpliceLocLast c) = c
lastContig (SpliceLocPrev _ n) = lastContig n

reloffsets :: SpliceLoc -> [Pos.Offset]
reloffsets = init . scanl (\off c -> off + CLoc.length c) 0 . contigs

-- | Returns the length of the region
length :: SpliceLoc -> Pos.Offset
length = foldl' (\len c -> len + CLoc.length c) 0 . contigs

-- | The bounds of a sequence location.  This is a pair consisting of
-- the lowest and highest sequence offsets covered by the region.  The
-- bounds ignore the strand of the sequence location, and the first
-- element of the pair will always be lower than the second.  Even if
-- the positions in the location do not run monotonically through the
-- location, the overall lowest and highest sequence offsets are returned.
bounds :: SpliceLoc -> (Pos.Offset, Pos.Offset)
bounds = (minimum *** maximum) . unzip . map CLoc.bounds . contigs

-- | Sequence position of the start of the location.  This is the 5'
-- end on the location strand, which will have a higher offset than
-- 'endPos' if the location is on the 'RevCompl' strand.
startPos :: SpliceLoc -> Pos.Pos
startPos = CLoc.startPos . firstContig

-- | Sequence position of the end of the location, as described in 'startPos'.
endPos :: SpliceLoc -> Pos.Pos
endPos = CLoc.endPos . lastContig

{-#


-- | Extract the nucleotide 'SeqData.SeqData' for the sequence location.  If
-- any part of the location lies outside the bounds of the sequence,
-- an error results.
seqDataM :: (MonadError m, Error (ErrorType m), SeqData.SeqData s, Stranded s) => s -> Loc -> m s
seqDataM sequ (Loc contigs)
    = liftM SeqData.concat $ mapM (CLoc.seqDataM sequ) contigs

-- | As 'seqData', extract the nucleotide subsequence for the
-- location.  Any positions in the location lying outside the bounds
-- of the sequence are returned as @N@ rather than producing an error.
seqDataPadded :: (SeqData.SeqData s, Stranded s) => s -> Loc -> s
seqDataPadded sequ (Loc contigs)
    = SeqData.concat $ map (CLoc.seqDataPadded sequ) contigs

-- | Given a sequence position and a sequence location relative to the
-- same sequence, compute a new position representing the original
-- position relative to the subsequence defined by the location.  If
-- the sequence position lies outside of the sequence location,
-- @Nothing@ is returned; thus, the offset of the new position will
-- always be in the range @[0, length cloc - 1]@.
--
-- When the sequence positions in the location are not monotonic,
-- there may be multiple possible posInto solutions.  That is, if the
-- same outer sequence position is covered by two different contiguous
-- blocks of the location, then it would have two possible sequence
-- positions relative to the location. In this case, the position
-- 5\'-most in the location orientation is returned.
posInto :: Pos.Pos -> Loc -> Maybe Pos.Pos
posInto seqpos (Loc contigs) = posIntoContigs seqpos contigs

posIntoContigs :: Pos.Pos -> [CLoc.ContigLoc] -> Maybe Pos.Pos
posIntoContigs _ [] = Nothing
posIntoContigs seqpos (contig@(CLoc.ContigLoc _ len _):rest)
    = case CLoc.posInto seqpos contig of
        just@(Just _) -> just
        Nothing -> liftM (flip Pos.slide len) $ posIntoContigs seqpos rest

-- | Given a sequence location and a sequence position within that
-- location, compute a new position representing the original position
-- relative to the outer sequence.  If the sequence position lies
-- outside the location, @Nothing@ is returned.  
-- 
-- This function inverts 'posInto' when the sequence position lies
-- within the position is actually within the location.  Due to the
-- possibility of redundant location-relative positions for a given
-- absolute position, 'posInto' does not necessary invert 'posOutof'
posOutof :: Pos.Pos -> Loc -> Maybe Pos.Pos
posOutof pos (Loc contigs) = posOutofContigs pos contigs

posOutofContigs :: Pos.Pos -> [CLoc.ContigLoc] -> Maybe Pos.Pos
posOutofContigs _ [] = Nothing
posOutofContigs seqpos (contig@(CLoc.ContigLoc _ len _):rest)
    = case CLoc.posOutof seqpos contig of
        just@(Just _) -> just
        Nothing -> posOutofContigs (Pos.slide seqpos $ negate len) rest

-- | Given a general disjoint sequence location and a contiguous
-- location on the same coordinate system, compute a new contiguous
-- location representing the original one relative to the coordinates
-- of the enclosing disjoint location.  The contiguous location must
-- lie entirely within one contiguous component of the enclosing
-- disjoint location.

clocInto :: CLoc.ContigLoc -> Loc -> Maybe CLoc.ContigLoc
clocInto subcloc = listToMaybe . catMaybes . map into . clocsAndOffsets
    where into (cloc, off) = liftM (CLoc.slide off) . CLoc.clocInto subcloc $ cloc

clocOutof :: CLoc.ContigLoc -> Loc -> Maybe Loc
clocOutof (CLoc.ContigLoc pos5 len cstrand) (Loc clocs)
    = liftM (stranded cstrand . Loc) . regionOutof pos5 len $ clocs

regionOutof :: Int64 -> Int64 -> [CLoc.ContigLoc] -> Maybe [CLoc.ContigLoc]
regionOutof _ _ [] = Nothing
regionOutof pos5 len (cloc0:rest)
    | pos5 < 0 = Nothing
    | len < 0 = error $ "regionOutof: input, " ++ show (pos5, len, cloc0)
    | pos5 >= len0 = regionOutof (pos5 - len0) len rest
    | (pos5 + len <= len0) = case CLoc.clocOutof (CLoc.ContigLoc pos5 len Fwd) cloc0 of
                               Just out0 -> Just [out0]
                               Nothing -> error $ "regionOutof: final bounds failure, " ++ show (pos5, len, cloc0)
    | otherwise = let outlen0 = len0 - pos5
                  in case CLoc.clocOutof (CLoc.ContigLoc pos5 outlen0 Fwd) cloc0 of
                       Just out0 -> liftM (out0 :) . regionOutof 0 (len - outlen0) $ rest
                       Nothing -> error $ "regionOutof: internal bounds failure, " ++ show (pos5, len, outlen0, cloc0)
    where len0 = CLoc.length cloc0

-- | Returns a sequence location produced by extending the original
-- location on each end, based on a pair of (/5\' extension/, /3\'
-- extension/).  These add contiguous positions to the 5\' and 3\'
-- ends of the original location.  The 5\' extension is applied to the
-- 5\' end of the location on the location strand; if the location is
-- on the 'RevCompl' strand, the 5\' end will have a higher offset
-- than the 3\' end and this offset will increase by the amount of the
-- 5\' extension.  Similarly, the 3\' extension is applied to the 3\'
-- end of the location.
extend :: (Int64, Int64) -- ^ (5' extension, 3' extension)
       -> Loc -> Loc
extend _ (Loc []) = error "extendLoc on zero-contig Loc"
extend (ext5, ext3) (Loc contigs) = Loc $ case extendContigs3 contigs of
                                            [] -> error "Empty contigs after extendContigs3"
                                            (cfirst:crest) -> (CLoc.extend (ext5, 0) cfirst):crest
    where extendContigs3 [] = error "Empty contigs in extendContigs3"
          extendContigs3 [clast] = [CLoc.extend (0, ext3) clast]
          extendContigs3 (contig:crest) = contig : extendContigs3 crest

-- | Returns @True@ when  a sequence position lies within a sequence
-- location on the same sequence, and occupies the same strand.
isWithin :: Pos.Pos -> Loc -> Bool
isWithin seqpos (Loc contigs) = or $ map (CLoc.isWithin seqpos) contigs

overlappingContigs :: Loc -> Loc -> [(CLoc.ContigLoc, CLoc.ContigLoc)]
overlappingContigs (Loc contigs1) (Loc contigs2) 
    = filter (uncurry CLoc.overlaps) [(c1, c2) | c1 <- contigs1, c2 <- contigs2 ]

-- | Returns @True@ when two sequence locations overlap at any
-- position.
overlaps :: Loc -> Loc -> Bool
overlaps l1 l2 = not $ null $ overlappingContigs l1 l2

-- | Display a human-friendly, zero-based representation of a sequence location.
displayLoc :: Loc -> BS.ByteString
displayLoc (Loc contigs) = BS.intercalate (BS.singleton ';') $ map display contigs

-- | Little 'Parsing.ParseBS' parser for 'display' format
-- representations of a 'Loc' location
parser :: Parsing.ParseBS Loc
parser = Loc <$> Parsing.mapSplit undisplayParser ';'
 #-}