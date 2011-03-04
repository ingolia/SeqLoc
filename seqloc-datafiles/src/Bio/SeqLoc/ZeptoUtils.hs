{-# LANGUAGE OverloadedStrings #-}

module Bio.SeqLoc.ZeptoUtils
       where

import Control.Applicative
import qualified Data.ByteString as BSW
import qualified Data.ByteString.Char8 as BS
import Data.ByteString.Internal (c2w)

import qualified Data.Attoparsec.Char8 as AP (isDigit_w8)
import qualified Data.Attoparsec.Zepto as ZP

import Debug.Trace

import Bio.SeqLoc.Strand

strand :: ZP.Parser Strand
strand = ((ZP.string "+\t" *> return Fwd) <|>
          (ZP.string "-\t" *> return RevCompl))

decfield :: (Integral a) => ZP.Parser a
decfield = decimal <* ZP.string "\t"

field :: ZP.Parser BS.ByteString
field =  ZP.takeWhile (/= c2w '\t') <* ZP.string "\t"

dropField :: ZP.Parser ()
dropField = ZP.takeWhile (/= c2w '\t') *> ZP.string "\t" *> return ()

decimal :: (Integral a) => ZP.Parser a
decimal = decode <$> ZP.takeWhile (AP.isDigit_w8)
  where decode = fromIntegral . BSW.foldl' step (0 :: Int)
        step a w = a * 10 + fromIntegral (w - 48)
