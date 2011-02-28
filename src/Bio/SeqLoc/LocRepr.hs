module Bio.SeqLoc.LocRepr
       ( LocRepr(..)
       , unreprMaybe, unreprEither, unreprErr
       )       
       where

import qualified Data.ByteString.Char8 as BS

import qualified Data.Attoparsec as AP

class LocRepr l where
  repr :: l -> BS.ByteString
  unrepr :: AP.Parser l
  
unreprMaybe :: (LocRepr l) => BS.ByteString -> Maybe l
unreprMaybe = AP.maybeResult . flip AP.feed BS.empty . AP.parse unrepr

unreprEither :: (LocRepr l) => BS.ByteString -> Either String l
unreprEither = AP.eitherResult . flip AP.feed BS.empty . AP.parse unrepr

unreprErr :: (LocRepr l) => BS.ByteString -> l
unreprErr = either error id . AP.eitherResult . flip AP.feed BS.empty . AP.parse unrepr