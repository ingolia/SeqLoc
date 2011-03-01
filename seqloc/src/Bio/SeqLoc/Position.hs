{-# LANGUAGE FlexibleContexts, GeneralizedNewtypeDeriving, OverloadedStrings #-}

{-| Data type for a sequence position.

Zero-based 'Offset' / 'Int64' indices are used throughout, to
facilitate direct use of indexing functions on 'SeqData'.
-}

module Bio.SeqLoc.Position ( 
  -- * Sequence positions
  Offset(..)
  , Pos(..)

  -- * Manipulating positions
  , slide

  -- * Extracting sequences
  , atOffset, atPos                           
  )
    where 

import Control.Applicative
import Control.Monad (liftM)
import qualified Data.ByteString as BSW
import qualified Data.ByteString.Char8 as BS
import qualified Data.ListLike as LL
import Data.Word (Word8)

import qualified Data.Attoparsec.Char8 as AP (isDigit_w8)
import qualified Data.Attoparsec.Zepto as ZP

import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.Strand

-- | Unstranded offset in a sequence
newtype Offset = Offset { unOffset :: Int } deriving (Eq, Ord, Show, Read, Num, Real, Enum, Integral)

instance LocRepr Offset where
  repr = BS.pack . show . unOffset
  unrepr = (negate <$> (ZP.string "-" *> decimal)) <|> (ZP.string "+" *> decimal) <|> decimal
    where decimal = Offset . BSW.foldl' step 0 <$> ZP.takeWhile AP.isDigit_w8
          step a w = a * 10 + fromIntegral (w - 48)

-- | Stranded position in a sequence
data Pos = Pos { offset :: !Offset -- ^ 0-based index of the position
               , strand :: !Strand -- ^ Strand of the position
               }
              deriving (Eq, Ord, Show, Read)
  
instance Stranded Pos where
  revCompl (Pos off str) = Pos off (revCompl str)
  
instance LocRepr Pos where
  repr (Pos off str) = BS.concat [ repr off, repr str ]
  unrepr = Pos <$> unrepr <*> unrepr

-- | Returns a position resulting from sliding the original position
-- along the sequence by a specified offset.  A positive offset will
-- move the position away from the 5\' end of the forward stand of the
-- sequence regardless of the strand of the position itself.  Thus,
-- 
-- > slide (revCompl pos) off == revCompl (slide pos off)
slide :: Pos -> Offset -> Pos
slide (Pos off str) doff = Pos (off + doff) str

-- | Extract 'Just' the element at a specific sequence offset, or
-- 'Nothing' if the offset is outside the bounds of the the region.
atOffset :: (LL.ListLike s c) => s -> Offset -> Maybe c
atOffset s (Offset off) | off < 0 = Nothing
                        | otherwise = case LL.drop off s of
                          trimmed | LL.null trimmed -> Nothing
                                  | otherwise -> Just . LL.head $ trimmed


-- | Extract 'Just' the item at a specific sequence position, or
-- 'Nothing' if the position lies outside the bounds of the sequence.
atPos :: (LL.ListLike s c, Stranded c) => s -> Pos -> Maybe c
atPos sequ (Pos off str) = liftM (stranded str) . atOffset sequ $ off

{-# SPECIALIZE atPos :: String -> Pos -> Maybe Char #-}
{-# SPECIALIZE atPos :: BS.ByteString -> Pos -> Maybe Word8 #-}