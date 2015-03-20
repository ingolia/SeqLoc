{-# LANGUAGE GeneralizedNewtypeDeriving, TypeSynonymInstances, FlexibleInstances #-}
{-| Efficient lookup based on potential location overlap

Collection of objects with locations, potentially on named sequences,
that can be queried to recover all objects whose location could
potentially overlap a query location.

|-}

module Bio.SeqLoc.LocMap (
  -- * Mapping objects to locations
  LocMap
  , emptyLM, insertLoc
  , queryLoc

  -- * Mapping objects to locations on named sequneces
  , SeqLocMap
  , emptySLM, insertSeqLoc
  , querySeqLoc

    -- * Mapping transcripts to their locations
  , transcriptSeqLocMap
    
  -- * Generalization of objects with locations
  , Locatable(..)
  , WithLocation(..)
  , locatableSeqLocMap
  , queryLocatable, queryLocCompatible
    
  )
       where

import qualified Data.HashMap.Strict as HM
import Data.List
-- import Data.Maybe

import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import qualified Bio.SeqLoc.ShiftedVector as ShV

-- | Mapping objects to sequence locations in the sense of `Loc.Location`.
data LocMap a = LocMap { binSize :: !Pos.Offset,
                         bins :: !(ShV.ShiftedVector [a]) } deriving (Show)

-- | Mapping objects to locations on named sequences 
data SeqLocMap a = SeqLocMap { slmBinSize :: !Pos.Offset, locmaps :: !(HM.HashMap SeqLabel (LocMap a)) } deriving (Show)

-- | Object with a genomic location expressed by `SpliceSeqLoc`.
class Locatable o where
  locate :: o -> SpliceSeqLoc

instance Locatable ContigSeqLoc where
  locate (OnSeq ref loc) = OnSeq ref (SpLoc.singleton loc)

instance Locatable SpliceSeqLoc where
  locate = id

instance Locatable Transcript where
  locate = location

-- | Simple representation of a `Locatable` object as an arbitrary
-- object adjoined with a location.
data WithLocation a = WithLocation { withoutLocation :: !a, withLocate :: !SpliceSeqLoc } deriving (Eq)

instance Locatable (WithLocation a) where
  locate = withLocate

-- | Create an empty object / location map.
--
-- Specify a characteristic size for efficient queries, which depends
-- on the underlying genome. Smaller sizes yield fewer false
-- candidates but greater memory usage and potentially slower queries.
emptyLM :: Pos.Offset -> LocMap a
emptyLM bsz = LocMap { binSize = bsz, bins = ShV.empty }

-- | Insert a new object / location pair and reutrn the updated map
insertLoc :: (Loc.Location l) => l -> a -> LocMap a -> LocMap a
insertLoc l x lm0 = lm0 { bins = ShV.modifySome (bins lm0) (locBins l lm0) (x :)  }

-- | Retrieve a list of objects that could potentially overlap the
-- query location.
--
-- Some objects may not actually overlap the query location. Some
-- objects may appear more than once in the list. However, no object
-- whose location truly overlaps the query will be missing from the
-- list.
--
-- Overlap is defined by one or more nucleotides in common between the
-- bounds of the two locations, regardless of the relative strands of
-- the query and the object location.
queryLoc :: (Loc.Location l) => l -> LocMap a -> [a]
queryLoc l lm = concat [ (bins lm) ShV.!? b | b <- locBins l lm ]

locBins :: (Loc.Location l) => l -> LocMap a -> [Int]
locBins l lm = let (start, end) = Loc.bounds l
                   binlow = fromIntegral $ start `div` binSize lm
                   binhigh = fromIntegral $ end `div` binSize lm
               in [binlow..binhigh]

-- | Create an empty object / location map
--
-- Specify a characteristic size as per `emptyLM`.
emptySLM :: Pos.Offset -> SeqLocMap a
emptySLM bsz = SeqLocMap { slmBinSize = bsz, locmaps = HM.empty }

-- | Insert a new object / location pair and reutrn the updated map
insertSeqLoc :: (Loc.Location l) => OnSeq l -> a -> SeqLocMap a -> SeqLocMap a
insertSeqLoc sl x slm0 = let lm0 = HM.lookupDefault (emptyLM $ slmBinSize slm0) (onSeqLabel sl) (locmaps slm0)
                             lm' = insertLoc (unOnSeq sl) x lm0
                         in slm0 { locmaps = HM.insert (onSeqLabel sl) lm' (locmaps slm0) }

-- | Retrieve a list of objects that could potentially overlap the
-- query location.
--
-- Some objects may not actually overlap the query location. Some
-- objects may appear more than once in the list. However, no object
-- whose location truly overlaps the query will be missing from the
-- list.
--
-- Overlap is defined by one or more nucleotides in common between the
-- bounds of the two locations, regardless of the relative strands of
-- the query and the object location, as well as the same name for the
-- underlying reference sequence.
--
-- To retrieve objects and test for actual overlap see the `Locatable`
-- interface and `queryLocatable` or `queryLocCompatible`.
querySeqLoc :: (Loc.Location l) => OnSeq l -> SeqLocMap a -> [a]
querySeqLoc sl slm = maybe [] (queryLoc (unOnSeq sl)) $
                     HM.lookup (onSeqLabel sl) (locmaps slm)

-- | Construct a mapping from transcripts to their genomic locations.
transcriptSeqLocMap :: Pos.Offset -> [Transcript] -> SeqLocMap Transcript
transcriptSeqLocMap bsz = foldl' insertTrx (emptySLM bsz)
  where insertTrx slm0 t = insertSeqLoc (location t) t slm0

-- | Construct a mapping from a general `Locatable` object to its
-- genomic location
locatableSeqLocMap :: (Locatable l) => Pos.Offset -> [l] -> SeqLocMap l
locatableSeqLocMap bsz = foldl' insertLocable (emptySLM bsz)
  where insertLocable slm0 t = insertSeqLoc (locate t) t slm0

-- | Retrieve a list of `Locatable` objects whose overall, contiguous
-- genomic coordinates intersect at any position the genomic interval
-- spanned by the specified `Loc.Location`. This does not require that
-- the spliced structure of the query is a subset of the spliced
-- structure of the `Locatable` nor that the query location lie
-- entirely within the hit location (contrast with
-- `queryLocCompatible`).
--
-- When a strand argument is given, restrict matches to those lying on
-- the same strand as the query location, for `Just Plus`, or the
-- opposite strand, for `Just Minus`.
queryLocatable :: (Locatable o, Loc.Location l) => Maybe Strand -> OnSeq l -> SeqLocMap o -> [o]
queryLocatable mstr qyloc = filter realHit . querySeqLoc qyloc
  where realHit sb = let sbloc = locate sb
                         (qylb, qyub) = Loc.bounds . unOnSeq $ qyloc
                         (sblb, sbub) = Loc.bounds . unOnSeq $ sbloc
                         qystr = Loc.strand . unOnSeq $ qyloc
                         sbstr = Loc.strand . unOnSeq $ sbloc
                     in (qylb <= sbub) && (qyub >= sblb) && (qystr `strandCompat` sbstr)
        strandCompat = case mstr of
          Nothing -> \_ _ -> True
          Just Plus -> (==)
          Just Minus -> (/=)

-- | Retrieve a list of `Locatable` objects whose spliced structure
-- contains the query location specifically.
queryLocCompatible :: (Locatable o) => Maybe Strand -> SpliceSeqLoc -> SeqLocMap o -> [o]
queryLocCompatible mstr qyloc = filter compatHit . querySeqLoc qyloc
  where compatHit sb = maybe False strandCompat $ unOnSeq qyloc `SpLoc.locInto` (unOnSeq . locate $ sb)
        strandCompat = case mstr of
          Nothing -> const True
          Just s -> (== s) . Loc.strand 
