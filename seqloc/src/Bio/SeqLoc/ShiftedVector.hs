{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns #-}
module Bio.SeqLoc.ShiftedVector
       where

import Prelude hiding (length)
import Control.Arrow
import Data.Maybe
import Data.Monoid
import qualified Data.Vector as V

data ShiftedVector a = ShiftedVector { zerois :: !Int, nullis :: !a, vector :: !(V.Vector a) } deriving (Show)

empty :: (Monoid a) => ShiftedVector a
empty = ShiftedVector { zerois = 0, nullis = mempty, vector = V.empty }

emptyZ :: a -> ShiftedVector a
emptyZ z = ShiftedVector { zerois = 0, nullis = z, vector = V.empty }

singleton :: (Monoid a) => Int -> a -> ShiftedVector a
singleton i x = ShiftedVector { zerois = i, nullis = mempty, vector = V.singleton x }

replicate :: (Monoid a) => Int -> Int -> a -> ShiftedVector a
replicate i0 n x = ShiftedVector { zerois = i0, nullis = mempty, vector = V.replicate n x }

length :: ShiftedVector a -> Int
length = V.length . vector

null :: ShiftedVector a -> Bool
null = V.null . vector

start :: ShiftedVector a -> Int
start = zerois

end :: ShiftedVector a -> Int
end sv = zerois sv + length sv - 1

(!?) :: ShiftedVector a -> Int -> a
(!?) sv i = fromMaybe (nullis sv) $ (vector sv) V.!? (i - zerois sv)

(//) :: ShiftedVector a -> [(Int, a)] -> ShiftedVector a
(//) sv0 ixs | Prelude.null ixs = sv0
             | otherwise = updateUnsafe sv'
  where ilow = minimum . map fst $ ixs
        ihigh = maximum . map fst $ ixs
        sv' = ensureLow ilow . ensureHigh ihigh $ sv0
        updateUnsafe sv = sv { vector = vector sv V.// jxs }
          where jxs = map (first $ subtract (zerois sv)) ixs        

modifySome :: ShiftedVector a -> [Int] -> (a -> a) -> ShiftedVector a
modifySome sv0 is f | Prelude.null is = sv0
                    | otherwise = modifyUnsafe sv'
  where ilow = minimum is
        ihigh = maximum is
        sv' = ensureLow ilow . ensureHigh ihigh $ sv0
        modifyUnsafe sv = let js = map (subtract (zerois sv)) is
                              ys = [ f (sv !? j) | j <- js ]
                          in sv { vector = vector sv V.// (zip js ys) }

ensureLow :: Int -> ShiftedVector a -> ShiftedVector a
ensureLow lb sv0 = case zerois sv0 - lb of
  down | down <= 0 -> sv0
       | otherwise ->  let !d = max down ((V.length . vector $ sv0) `div` 2)
                       in sv0 { zerois = zerois sv0 - d, vector = (V.replicate d $ nullis sv0) V.++ vector sv0 }

ensureHigh :: Int -> ShiftedVector a -> ShiftedVector a
ensureHigh ub sv0 = case ub - end sv0 of
  up | up <= 0 -> sv0
     | otherwise -> let !u = max up ((V.length . vector $ sv0) `div` 2)
                    in sv0 { vector = vector sv0 V.++ (V.replicate u $ nullis sv0) }

