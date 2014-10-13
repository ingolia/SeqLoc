{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.SeqLoc.LocMap
       where

import qualified Data.HashMap.Strict as HM
import qualified Data.Vector as V

import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos

import qualified Bio.SeqLoc.ShiftedVector as ShV

data LocMap a = LocMap { binSize :: !Pos.Offset,
                         bins :: !(ShV.ShiftedVector [a]) } deriving (Show)

data SeqLocMap a = SeqLocMap { slmBinSize :: !Pos.Offset, locmaps :: !(HM.HashMap SeqLabel (LocMap a)) } deriving (Show)

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

emptySLM :: Pos.Offset -> SeqLocMap a
emptySLM bsz = SeqLocMap { slmBinSize = bsz, locmaps = HM.empty }

insertSeqLoc :: (Loc.Location l) => OnSeq l -> a -> SeqLocMap a -> SeqLocMap a
insertSeqLoc sl x slm0 = let lm0 = HM.lookupDefault (emptyLM $ slmBinSize slm0) (onSeqLabel sl) (locmaps slm0)
                             lm' = insertLoc (unOnSeq sl) x lm0
                         in slm0 { locmaps = HM.insert (onSeqLabel sl) lm' (locmaps slm0) }

querySeqLoc :: (Loc.Location l) => OnSeq l -> SeqLocMap a -> [a]
querySeqLoc sl slm = maybe [] (queryLoc (unOnSeq sl)) $
                     HM.lookup (onSeqLabel sl) (locmaps slm)
