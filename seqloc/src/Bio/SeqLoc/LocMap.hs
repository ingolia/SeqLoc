{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.SeqLoc.LocMap
       where

import qualified Data.Vector as V

import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos

import qualified Bio.SeqLoc.ShiftedVector as ShV

data LocMap a = LocMap { binSize :: !Pos.Offset,
                         bins :: !(ShV.ShiftedVector [a]) } deriving (Show)

emptyLM :: Pos.Offset -> LocMap a
emptyLM bsz = LocMap { binSize = bsz, bins = ShV.empty }

insertLoc :: (Loc.Location l) => l -> a -> LocMap a -> LocMap a
insertLoc l x lm0 = lm0 { bins = ShV.modifySome (bins lm0) (locBins l lm0) (x :)  }

queryLoc :: (Loc.Location l) => l -> LocMap a -> [a]
queryLoc l lm = concat [ (bins lm) ShV.!? b | b <- locBins l lm ]

locBins :: (Loc.Location l) => l -> LocMap a -> [Int]
locBins l lm = let (start, end) = Loc.bounds l
                   binlow = fromIntegral $ start `div` binSize lm
                   binhigh = fromIntegral $ end `div` binSize lm
               in [binlow..binhigh]
