{-# LANGUAGE FlexibleContexts, GeneralizedNewtypeDeriving, OverloadedStrings #-}

{-| Data type for a sequence position.

Zero-based 'Offset'indices are used throughout, to facilitate direct
use of indexing functions on 'SeqData'.  

-}

module Bio.SeqLoc.Position ( 
  -- * Sequence positions
  Offset(..)
  , Pos(..)

  -- * Manipulating positions
  , slide

  -- * Extracting sequences
  , atPos                           
  )
    where 

import Control.Applicative
import Control.Arrow
import Control.Monad (liftM)
import qualified Data.ByteString.Char8 as BS

import Bio.Core.Sequence
import Bio.Core.Strand

import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.SeqLike
import Bio.SeqLoc.Strand

-- | Stranded position in a sequence
data Pos = Pos { offset :: !Offset -- ^ 0-based index of the position
               , strand :: !Strand -- ^ Strand of the position
               }
              deriving (Eq, Ord, Show)
  
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

-- | Extract 'Just' the item at a specific sequence position, or
-- 'Nothing' if the position lies outside the bounds of the sequence.
atPos :: (SeqLike s) => s -> Pos -> Maybe Char
atPos sequ (Pos off str) = liftM (stranded str) . ntAt sequ $ off

{-# SPECIALIZE atPos :: String -> Pos -> Maybe Char #-}
{-# SPECIALIZE atPos :: BS.ByteString -> Pos -> Maybe Char #-}
