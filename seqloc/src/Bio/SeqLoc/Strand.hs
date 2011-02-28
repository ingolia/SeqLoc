{-# LANGUAGE OverloadedStrings, TypeSynonymInstances #-}

{-| Utilities for manipulating nucleotide sequences and locations on
nucleotide sequences that occur on a forward or a reverse-complement
strand.

-}

module Bio.SeqLoc.Strand ( Strand(..)
                         , compl
                         , Stranded(..), stranded
                         )
    where 

import Control.Applicative
import Data.ByteString.Internal (c2w, w2c)
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import Data.Word (Word8)

import qualified Data.Attoparsec.Zepto as ZP

import Bio.SeqLoc.LocRepr

-- | Complement of a nucleotide character, swap A/T and G/C preserving
-- case and leave all other characters unchanged.
compl :: Char -> Char
compl 'a' = 't'
compl 'c' = 'g'
compl 'g' = 'c'
compl 't' = 'a'
compl 'A' = 'T'
compl 'C' = 'G'
compl 'G' = 'C'
compl 'T' = 'A'
compl ch  = ch

-- | Sequence strand
data Strand = Fwd | RevCompl deriving (Eq, Ord, Show, Read, Bounded, Enum)

instance LocRepr Strand where
  repr Fwd = BS.pack "(+)"
  repr RevCompl = BS.pack "(-)"
  unrepr = (ZP.string "(+)" *> return Fwd) <|>
           (ZP.string "(-)" *> return RevCompl)
                             
-- | A nucleotide sequence or location on a nucleotide sequence that
--   lies on a specific strand and has an orientation.
class Stranded s where
    revCompl :: s -> s

-- | Convert the orientation of a 'Stranded' thing based on a
--   specified 'Strand'
stranded :: (Stranded s) => Strand -> s -> s
stranded Fwd      = id
stranded RevCompl = revCompl

instance Stranded Strand where
    revCompl Fwd      = RevCompl
    revCompl RevCompl = Fwd

instance Stranded Char where
    revCompl = compl

instance Stranded Word8 where
  revCompl = c2w . compl . w2c

instance Stranded String where
    revCompl = reverse . map compl

instance Stranded LBS.ByteString where
    revCompl = LBS.reverse . LBS.map compl

instance Stranded BS.ByteString where
    revCompl = BS.reverse . BS.map compl
    