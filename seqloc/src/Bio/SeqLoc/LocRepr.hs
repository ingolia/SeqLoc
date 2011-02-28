module Bio.SeqLoc.LocRepr
       ( LocRepr(..)
       , unreprMaybe, unreprEither, unreprErr
       )       
       where

import qualified Data.ByteString.Char8 as BS

import qualified Data.Attoparsec.Zepto as ZP

class LocRepr l where
  repr :: l -> BS.ByteString
  unrepr :: ZP.Parser l
  
unreprMaybe :: (LocRepr l) => BS.ByteString -> Maybe l
unreprMaybe = either (const Nothing) Just . ZP.parse unrepr

unreprEither :: (LocRepr l) => BS.ByteString -> Either String l
unreprEither = ZP.parse unrepr

unreprErr :: (LocRepr l) => BS.ByteString -> l
unreprErr = either error id . ZP.parse unrepr