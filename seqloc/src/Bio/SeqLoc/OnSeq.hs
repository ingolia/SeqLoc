{-# LANGUAGE FlexibleInstances, TypeFamilies, FlexibleContexts, GeneralizedNewtypeDeriving #-}

{-| 

Data types for sequence locations and sequence positions associated
with specific, named sequences.

-}

module Bio.SeqLoc.OnSeq ( 
  SeqLabel(..), toSeqLabel, unSeqLabel
  
  , OnSeq(..)
  
  -- * Positions on named sequences
  , SeqOffset, SeqPos
  
  -- * Contiguous location spans on named sequences
  , ContigSeqLoc

    -- * Arbitrary location spans on named sequences
  , SpliceSeqLoc
    
  , andSameSeq
    
  )
    where 

import Control.Applicative
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import Data.ByteString.Internal (c2w)

import Data.Hashable

import qualified Data.Attoparsec.Zepto as ZP

import Bio.Core.Sequence

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand

toSeqLabel :: BS.ByteString -> SeqLabel
toSeqLabel = SeqLabel . LBS.fromChunks . (: [])

unSeqLabel :: SeqLabel -> BS.ByteString
unSeqLabel = BS.concat . LBS.toChunks . unSL

data OnSeq s = OnSeq { onSeqLabel :: !SeqLabel, unOnSeq :: !s } deriving (Eq, Ord)

at :: BS.ByteString
at = BS.singleton '@'

instance Hashable SeqLabel where
  hashWithSalt salt = hashWithSalt salt . unSL
  hash = hash . unSL

instance Stranded s => Stranded (OnSeq s) where
  revCompl (OnSeq name obj) = OnSeq name (revCompl obj)

instance LocRepr s => LocRepr (OnSeq s) where
  repr (OnSeq name obj) = BS.concat [ unSeqLabel name, at, repr obj ]
  unrepr = OnSeq <$> (toSeqLabel <$> ZP.takeWhile (/= c2w '@')) <* ZP.string at <*> unrepr

type SeqOffset = OnSeq Pos.Offset

-- | A position on a named sequence
type SeqPos = OnSeq Pos.Pos

-- | A location consisting of a contiguous span of positions on a
-- named sequence.
type ContigSeqLoc = OnSeq Loc.ContigLoc

-- | A general location, consisting of spans of sequence positions on
-- a specific, named sequence.
type SpliceSeqLoc = OnSeq SpLoc.SpliceLoc

andSameSeq :: (a -> b -> Bool) -> OnSeq a -> OnSeq b -> Bool
andSameSeq p (OnSeq n1 x) (OnSeq n2 y) | n1 == n2 = p x y
                                       | otherwise = False

instance BioSeq (OnSeq BS.ByteString) where
  seqlabel = onSeqLabel
  seqdata = SeqData . LBS.fromChunks . (: []) . unOnSeq
  seqlength = Offset . fromIntegral . BS.length . unOnSeq
