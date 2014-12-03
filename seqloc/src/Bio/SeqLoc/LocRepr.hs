{-# LANGUAGE OverloadedStrings #-}

module Bio.SeqLoc.LocRepr
       ( LocRepr(..), reprStr
       , unreprMaybe, unreprEither, unreprErr
       )       
       where

import Control.Applicative
import qualified Data.ByteString as BSW
import qualified Data.ByteString.Char8 as BS

import qualified Data.Attoparsec.Char8 as AP (isDigit_w8)
import qualified Data.Attoparsec.Zepto as ZP

import Bio.Core.Sequence
import Bio.Core.Strand

class LocRepr l where
  repr :: l -> BS.ByteString
  unrepr :: ZP.Parser l
  
reprStr :: (LocRepr l) => l -> String
reprStr = BS.unpack . repr
  
unreprMaybe :: (LocRepr l) => BS.ByteString -> Maybe l
unreprMaybe = either (const Nothing) Just . ZP.parse unrepr

unreprEither :: (LocRepr l) => BS.ByteString -> Either String l
unreprEither = ZP.parse unrepr

unreprErr :: (LocRepr l) => BS.ByteString -> l
unreprErr = either error id . ZP.parse unrepr

instance LocRepr Strand where
  repr Plus = BS.pack "(+)"
  repr Minus = BS.pack "(-)"
  unrepr = (ZP.string "(+)" *> return Plus) <|>
           (ZP.string "(-)" *> return Minus)

instance LocRepr Offset where
  repr = BS.pack . show . unOff
  unrepr = (negate <$> (ZP.string "-" *> decimal)) <|> (ZP.string "+" *> decimal) <|> decimal
    where decimal = Offset . BSW.foldl' step 0 <$> ZP.takeWhile AP.isDigit_w8
          step a w = a * 10 + fromIntegral (w - 48)

