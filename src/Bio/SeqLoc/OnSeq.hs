{-# LANGUAGE FlexibleInstances, TypeFamilies, FlexibleContexts, GeneralizedNewtypeDeriving #-}

{-| 

Data types for sequence locations and sequence positions associated
with specific, named sequences.

-}

module Bio.SeqLoc.OnSeq ( 
  SeqName(..)
  
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
import Data.String

import qualified Data.Attoparsec.Char8 as AP

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc

newtype SeqName = SeqName { unSeqName :: BS.ByteString } deriving (Show, Read, Eq, Ord, IsString)

data OnSeq s = OnSeq { onSeqName :: !SeqName, unOnSeq :: !s } deriving (Show, Read, Eq, Ord)

at :: BS.ByteString
at = BS.singleton '@'

instance LocRepr s => LocRepr (OnSeq s) where
  repr (OnSeq name obj) = BS.concat [ unSeqName name, at, repr obj ]
  unrepr = OnSeq <$> (SeqName <$> AP.takeWhile1 (/= '@')) <* AP.char '@' <*> unrepr

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
