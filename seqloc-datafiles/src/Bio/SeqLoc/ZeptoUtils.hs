{-# LANGUAGE OverloadedStrings #-}

module Bio.SeqLoc.ZeptoUtils
       where

import Control.Applicative
import qualified Data.ByteString as BSW
import qualified Data.ByteString.Char8 as BS
import Data.ByteString.Internal (c2w, w2c)
import Data.Word

import qualified Data.Attoparsec.Char8 as AP (isDigit_w8)
import qualified Data.Attoparsec.Zepto as ZP

import Bio.SeqLoc.Strand

firstfield :: ZP.Parser BS.ByteString
firstfield = ZP.takeWhile (/= c2w '\t')

strand :: ZP.Parser Strand
strand = ((ZP.string "\t+" *> return Plus) <|>
          (ZP.string "\t-" *> return Minus))

decfield :: (Integral a) => ZP.Parser a
decfield = ZP.string "\t" *> decimal

field :: ZP.Parser BS.ByteString
field =  ZP.string "\t" *> ZP.takeWhile (/= c2w '\t')

dropField :: ZP.Parser ()
dropField = ZP.string "\t" *> ZP.takeWhile (/= c2w '\t') *> return ()

decimal :: (Integral a) => ZP.Parser a
decimal = decode <$> ZP.takeWhile (AP.isDigit_w8)
  where decode = fromIntegral . BSW.foldl' step (0 :: Int)
        step a w = a * 10 + fromIntegral (w - 48)

decode :: (Integral a) => BSW.ByteString -> Either String a
decode str = either Left (Right . fromIntegral) . BSW.foldl' step (Right (0 :: Int)) $! str
  where step :: Either String Int -> Word8 -> Either String Int
        step err@(Left _) _ch = err
        step (Right a) w | w >= 48 && w < 58 = Right $! a * 10 + fromIntegral (w - 48)
        step _ w = Left $ "Bad character " ++ show (w2c w) ++ " in number " ++ show str

unlessAtEnd :: ZP.Parser a -> ZP.Parser (Maybe a)
unlessAtEnd p = do isAtEnd <- ZP.atEnd
                   if isAtEnd then return Nothing else fmap Just p
