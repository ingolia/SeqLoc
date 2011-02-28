{-# LANGUAGE FlexibleInstances #-}

module Bio.SeqLoc.SeqData
       ( SeqData(..)
       )
       where 

import Prelude hiding (length, concat)

import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import qualified Data.List as L
import Data.Maybe

class SeqData s where
  -- | Length of sequence data
  length :: (Integral n) => s -> n
  
  -- | 'Just' the nucleotide at a specified sequence data offset,
  -- given in 0-based coordinates, or 'Nothing' if the offset is
  -- beyond the bounds of the data
  ntAt :: (Integral n) => s -> n -> Maybe Char
  
  -- | 'Just' the nucleotides in subsequence of the sequence data, or
  -- 'Nothing' if the region extends beyond the bounds of the
  -- sequence.
  subseq :: (Integral n, Integral m) 
            => n -- ^ Starting position in 0-based coordinates
            -> m -- ^ Length
            -> s -- ^ Sequence data
            -> Maybe s
            
  -- | Nucleotides in a subsequence of the sequence data, padded with
  -- @N@ when the region extends beyond the bounds of the sequence.
  subseqPad :: (Integral n, Integral m)
               => n -- ^ Starting position in 0-based coordinates
               -> m -- ^ Length
               -> s -- ^ Sequence data
               -> s
  
  concat :: [s] -> s

instance SeqData [Char] where
  length = L.genericLength
  ntAt l p | p < 0 = Nothing
           | otherwise = listToMaybe . drop (fromIntegral p) $ l
  subseq = listSubseq
  subseqPad = listSubseqPad
  concat = L.concat
  
listSubseq :: (Integral n, Integral m) => n -> m -> [Char] -> Maybe [Char]
listSubseq start len l | start < 0 = Nothing
                       | otherwise = case take (fromIntegral len) . drop (fromIntegral start) $ l of 
                         ss | L.genericLength ss == len -> Just ss
                            | otherwise -> Nothing

listSubseqPad :: (Integral n, Integral m) => n -> m -> [Char] -> [Char]
listSubseqPad istart ilen sequ
  | start + len <= 0 = replicate len 'N'
  | start < 0 = replicate (negate start) 'N' ++ takePadded (len + start) sequ
  | otherwise = takePadded len . drop start $ sequ
    where takePadded sublen subsequ = case take sublen subsequ of
            subsubsequ | length subsubsequ == sublen -> subsubsequ
                       | otherwise -> subsubsequ ++ replicate (sublen - length subsubsequ) 'N'
          start = fromIntegral istart
          len = fromIntegral ilen

instance SeqData BS.ByteString where      
  length = fromIntegral . BS.length
  ntAt sequ ipos | pos >= 0 && pos < BS.length sequ = Just $ BS.index sequ pos
                 | otherwise = Nothing
                  where pos = fromIntegral ipos
  subseq = bsSubseq
  subseqPad = bsSubseqPad
  concat = BS.concat
  
bsSubseq :: (Integral n, Integral m) => n -> m -> BS.ByteString -> Maybe BS.ByteString
bsSubseq istart ilen sequ
  | start < 0 || start + len > BS.length sequ = Nothing
  | otherwise = Just . BS.take len . BS.drop start $ sequ
    where start = fromIntegral istart
          len = fromIntegral ilen

bsSubseqPad :: (Integral n, Integral m) => n -> m -> BS.ByteString -> BS.ByteString                
bsSubseqPad istart ilen sequ
  | start + len <= 0        = BS.replicate len 'N'
  | start >= BS.length sequ = BS.replicate len 'N'
  | start < 0  = BS.replicate (negate start) 'N' `BS.append` takePadded (len + start) sequ
  | otherwise = takePadded len $ BS.drop start sequ
    where takePadded sublen subsequ
            | sublen <= BS.length subsequ = BS.take sublen subsequ
            | otherwise = subsequ `BS.append` BS.replicate (sublen - BS.length subsequ) 'N'
          start = fromIntegral istart
          len = fromIntegral ilen

instance SeqData LBS.ByteString where
    length = fromIntegral . LBS.length
    ntAt = lbsNtAt
    subseq = lbsSubseq
    subseqPad = lbsSubseqPad
    concat = LBS.concat
    
lbsNtAt :: (Integral n) => LBS.ByteString -> n -> Maybe Char
lbsNtAt sequ pos | pos < 0 = Nothing
                 | otherwise = case LBS.drop (fromIntegral pos) sequ of
                   trimmed | LBS.null trimmed -> Nothing
                           | otherwise -> Just . LBS.head $ trimmed

lbsSubseq :: (Integral n, Integral m) => n -> m -> LBS.ByteString -> Maybe LBS.ByteString
lbsSubseq istart ilen sequ
  | start < 0 = Nothing
  | otherwise = case LBS.take len . LBS.drop start $ sequ of
    subsequ | LBS.length subsequ < len -> Nothing
            | otherwise -> Just subsequ
    where start = fromIntegral istart
          len = fromIntegral ilen
    
lbsSubseqPad :: (Integral n, Integral m) => n -> m -> LBS.ByteString -> LBS.ByteString
lbsSubseqPad istart ilen sequ
  | start + len <= 0 = LBS.replicate len 'N'
  | start < 0  = LBS.replicate (negate start) 'N' `LBS.append` takePadded (len + start) sequ
  | otherwise = takePadded len $ LBS.drop start sequ
    where takePadded sublen subsequ = case LBS.take sublen subsequ of
            subsubsequ | LBS.length subsubsequ == sublen -> subsubsequ
                       | otherwise -> subsubsequ `LBS.append` LBS.replicate (sublen - LBS.length subsubsequ) 'N'
          start = fromIntegral istart
          len = fromIntegral ilen