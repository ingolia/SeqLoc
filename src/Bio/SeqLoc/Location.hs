{-# LANGUAGE TypeFamilies, FlexibleContexts #-}
{-| Data type for a sequence location consiting of a contiguous range
of positions on the sequence.

Throughout, /sequence position/ refers to a 'Pos.Pos' which includes a
strand, as opposed to an /offset/, which refers to a 'Pos.Offset' with
no strand.

 -}

module Bio.SeqLoc.Location 
       ( -- * Sequence locations
         Location(..)
         -- * Contiguous sequence locations
       , ContigLoc(..), fromStartEnd, fromPosLen
         -- * Transforming locations
       , slide
         -- * Printing locations
       , displayContigLoc
       ) 
    where

import Prelude hiding (length)

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS

import qualified Data.Attoparsec.Char8 as AP

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import qualified Bio.SeqLoc.SeqData as SeqData

class Location l where
  length :: l -> Pos.Offset
  
  -- | The bounds of a sequence location.  This is a pair consisting
  -- of the lowest and highest sequence offsets covered by the region.
  -- The bounds ignore the strand of the sequence location, and the
  -- first element of the pair will always be lower than the second.
  bounds :: l -> (Pos.Offset, Pos.Offset)

  -- | Sequence position of the start of the location.  This is the 5'
  -- end on the location strand, which will have a higher offset than
  -- 'endPos' if the location is on the 'RevCompl' strand.
  startPos :: l -> Pos.Pos

  -- | Sequence position of the end of the location, as described in
  -- 'startPos'.
  endPos :: l -> Pos.Pos

  -- | Extract 'Just' the nucleotide 'SeqData' for the sequence
  -- location, or 'Nothing' if f any part of the location lies outside
  -- the bounds of the sequence.
  seqData :: (SeqData.SeqData s, Stranded s) => s -> l -> Maybe s
  
  -- | As 'seqData', extract the nucleotide subsequence for the
  -- location, but any positions in the location lying outside the
  -- bounds of the sequence are returned as @N@.
  seqDataPad :: (SeqData.SeqData s, Stranded s) => s -> l -> s
  
  -- | Given a sequence position and a sequence location relative to
  -- the same sequence, compute a new position representing the
  -- original position relative to the subsequence defined by the
  -- location.  If the sequence position lies outside of the sequence
  -- location, @Nothing@ is returned; thus, the offset of the new
  -- position will always be in the range @[0, length l - 1]@.
  posInto :: Pos.Pos -> l -> Maybe Pos.Pos
  
  -- | Given a sequence location and a sequence position within that
  -- location, compute a new position representing the original
  -- position relative to the outer sequence.  If the sequence
  -- position lies outside the location, @Nothing@ is returned.
  -- 
  -- This function inverts 'posInto' when the sequence position lies
  -- within the position is actually within the location.
  posOutof :: Pos.Pos -> l -> Maybe Pos.Pos  
  
  -- | For an enclosing location and a sublocation in the same
  -- coordinate system, find the image of the sublocation relative to
  -- the enclosing location.  For example, if the enclosing location
  -- spans (100, 150) and the sublocation is (110, 120) then
  -- 'clocInto' will be (10, 20).
  clocInto :: ContigLoc -> l -> Maybe ContigLoc
    
  -- | Returns a sequence location produced by finding the inverse
  -- image of a sublocation, with coordinates given relative to an
  -- enclosing location, in the coordinate system of the enclosing
  -- location.  For example, if the enclosing location spans (100,
  -- 150) and the sublocation is (10, 20) then 'clocOutof' will be
  -- (110, 120).
  clocOutof :: ContigLoc -> l -> Maybe l

  -- | Returns a sequence location produced by extending the original
  -- location on each end, based on a pair of (/5\' extension/, /3\'
  -- extension/).  The 5\' extension is applied to the 5\' end of the
  -- location on the location strand; if the location is on the
  -- 'RevCompl' strand, the 5\' end will have a higher offset than the
  -- 3\' end and this offset will increase by the amount of the 5\'
  -- extension.  Similarly, the 3\' extension is applied to the 3\'
  -- end of the location.
  extend :: (Pos.Offset, Pos.Offset) -> l -> l

  -- | Returns @True@ when a sequence offset lies within a sequence
  -- location on the same sequence
  offsetWithin :: Pos.Offset -> l -> Bool

  -- | Returns @True@ when a sequence position lies within a sequence
  -- location on the same sequence, and occupies the same strand.
  posWithin :: Pos.Pos -> l -> Bool

  -- | Returns @True@ when two sequence locations overlap at any
  -- position.
  contigOverlaps :: ContigLoc -> l -> Bool

  toContigs :: l -> [ContigLoc]

  overlaps :: (Location l1) => l -> l1 -> Bool
  overlaps l0 = any (\c1 -> contigOverlaps c1 l0) . toContigs

-- | Contiguous sequence location defined by a span of sequence
-- positions, lying on a specific strand of the sequence.
data ContigLoc = ContigLoc { offset5 :: !Pos.Offset   -- ^ The offset of the 5\' end of the location, as a 0-based index
                           , clocLength :: !Pos.Offset    -- ^ The length of the location
                           , strand :: !Strand    -- ^ The strand of the location
                           } deriving (Eq, Ord, Show)

instance Stranded ContigLoc where
  revCompl (ContigLoc seq5 len str) = ContigLoc seq5 len $ revCompl str

to :: BS.ByteString
to = BS.pack "to"

instance LocRepr ContigLoc where
  repr cloc = let (seq5, seq3) = bounds cloc 
              in BS.concat [ repr seq5, to, repr seq3, repr . strand $ cloc ]
  unrepr = fromBoundsStrand <$> unrepr <* AP.string to <*> unrepr <*> unrepr

instance Location ContigLoc where
  length = clocLength
  seqData sequ (ContigLoc seq5 len str) = liftM (stranded str) . (SeqData.subseq seq5 len) $ sequ
  seqDataPad sequ (ContigLoc seq5 len str) = (stranded str) . (SeqData.subseqPad seq5 len) $ sequ
  posInto = clocPosInto
  posOutof = clocPosOutof
  bounds (ContigLoc seq5 len _) = (seq5, seq5 + len - 1)
  startPos (ContigLoc seq5 len str) 
    = case str of
        Fwd      -> Pos.Pos seq5             str
        RevCompl -> Pos.Pos (seq5 + len - 1) str
  endPos (ContigLoc seq5 len str) 
    = case str of
        Fwd      -> Pos.Pos (seq5 + len - 1) str
        RevCompl -> Pos.Pos seq5             str
  clocInto = clocClocInto
  clocOutof = clocClocOutof
  extend = clocExtend
  offsetWithin off (ContigLoc seq5 len _)
    = (off >= seq5) && (off < seq5 + len)
  posWithin (Pos.Pos pos pStrand) (ContigLoc seq5 len cStrand) 
    = (pos >= seq5) && (pos < seq5 + len) && (cStrand == pStrand)
  contigOverlaps = clocOverlaps
  toContigs = (: [])

-- | Create a sequence location between 0-based starting and ending
-- bounds with a specified strand.
fromBoundsStrand :: Pos.Offset -> Pos.Offset -> Strand -> ContigLoc
fromBoundsStrand seq5 seq3 _ | seq3 < seq5 = error "Bio.SeqLoc.Location.fromBoundsStrand: seq3 < seq5"
fromBoundsStrand seq5 seq3 str = ContigLoc seq5 (1 + seq3 - seq5) str

-- | Create a sequence location lying between 0-based starting and
-- ending offsets.  When @start < end@, the location
-- be on the forward strand, otherwise it will be on the
-- reverse complement strand.
fromStartEnd :: Pos.Offset -> Pos.Offset -> ContigLoc
fromStartEnd start end
    | start < end = ContigLoc start (1 + end - start) Fwd
    | otherwise   = ContigLoc end   (1 + start - end) RevCompl

-- | Create a sequence location from the sequence position of the
-- start of the location and the length of the position.  The strand
-- of the location, and the direction it extends from the starting
-- position, are determined by the strand of the starting position.
fromPosLen :: Pos.Pos -> Pos.Offset -> ContigLoc
fromPosLen _                       len | len < 0 = error "Bio.SeqLoc.Location.fromPosLen: len < 0"
fromPosLen (Pos.Pos off5 Fwd)      len = ContigLoc off5               len Fwd
fromPosLen (Pos.Pos off3 RevCompl) len = ContigLoc (off3 - (len - 1)) len RevCompl

-- | Returns a location resulting from sliding the original location
-- along the sequence by a specified offset.  A positive offset will
-- move the location away from the 5\' end of the forward stand of the
-- sequence regardless of the strand of the location itself.  Thus,
-- 
-- > slide (revCompl cloc) off == revCompl (slide cloc off)
slide :: Pos.Offset -> ContigLoc -> ContigLoc
slide dpos (ContigLoc seq5 len str) = ContigLoc (seq5 + dpos) len str

clocPosInto :: Pos.Pos -> ContigLoc -> Maybe Pos.Pos
clocPosInto (Pos.Pos pos pStrand) (ContigLoc seq5 len cStrand)
    | pos < seq5 || pos >= seq5 + len = Nothing
    | otherwise = Just $ case cStrand of
                           Fwd      -> Pos.Pos (pos - seq5)              pStrand
                           RevCompl -> Pos.Pos (seq5 + len  - (pos + 1)) (revCompl pStrand)

clocPosOutof :: Pos.Pos -> ContigLoc -> Maybe Pos.Pos
clocPosOutof (Pos.Pos pos pStrand) (ContigLoc seq5 len cStrand)
    | pos < 0 || pos >= len = Nothing
    | otherwise = Just $ case cStrand of
                           Fwd      -> Pos.Pos (pos + seq5)             pStrand
                           RevCompl -> Pos.Pos (seq5 + len - (pos + 1)) (revCompl pStrand)

clocExtend :: (Pos.Offset, Pos.Offset) -> ContigLoc -> ContigLoc
clocExtend (ext5, ext3) (ContigLoc seq5 len str)
    = case str of
        Fwd -> ContigLoc (seq5 - ext5) (len + ext5 + ext3) str
        RevCompl -> ContigLoc (seq5 - ext3) (len + ext5 + ext3) str

clocClocInto :: ContigLoc -> ContigLoc -> Maybe ContigLoc
clocClocInto subcloc tocloc
    = case (posInto (startPos subcloc) tocloc, posInto (endPos subcloc) tocloc) of
        (Just start, Just _) -> Just (fromPosLen start (length subcloc))
        _ -> Nothing

clocClocOutof :: ContigLoc -> ContigLoc -> Maybe ContigLoc
clocClocOutof subcloc fromcloc
    = case (posOutof (startPos subcloc) fromcloc, posOutof (endPos subcloc) fromcloc) of
        (Just start, Just _) -> Just (fromPosLen start (length subcloc))
        _ -> Nothing

clocOverlaps :: ContigLoc -> ContigLoc -> Bool
clocOverlaps contig1 contig2
    = case (bounds contig1, bounds contig2) of
        ((low1, high1),(low2, high2)) -> (strand contig1 == strand contig2)
                                         && (low1 <= high2) && (low2 <= high1)

displayContigLoc :: ContigLoc -> BS.ByteString
displayContigLoc cloc = BS.concat [ BS.pack . show . Pos.offset . startPos $ cloc
                                  , BS.pack "to"
                                  , BS.pack . show . Pos.offset . endPos $ cloc
                                  ]

--parser :: Parsing.ParseBS ContigLoc
--parser = fromStartEnd <$>
--         (Parsing.integral <* Parsing.string "to") <*>
--         Parsing.integral
