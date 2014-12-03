{-# LANGUAGE BangPatterns, ExistentialQuantification #-}
module Main
       where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe

import Test.QuickCheck
import Test.QuickCheck.Monadic

import qualified Bio.SeqLoc.Bed as Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

main :: IO ()
main = mapM_ runTest tests

tests :: [Test]
tests = [ T "Write random bed"          test_writeRandomBed ]

test_writeRandomBed :: Property
test_writeRandomBed = monadicIO $
                      forAllM (listOf1 genTranscript) $ \trxs ->
                      run . BS.writeFile "tmpf.bed" . BS.unlines . map Bed.transcriptToBedStd $ trxs
                      

-- | 

genName :: Gen SeqLabel
genName = liftM (toSeqLabel . BS.pack) $ genNameLength >>= flip replicateM genNameChar
    where genNameLength = choose (1, 20)
          genNameChar = elements $ ['a'..'z'] ++ ['A'..'Z'] ++ ['0'..'9'] ++ "-_"

instance Arbitrary SeqLabel where
    arbitrary = genName
  
-- | Constrained position generators

genOffset :: Gen Pos.Offset
genOffset = do isneg <- arbitrary
               nnoff <- genNonNegOffset
               return $ (if isneg then negate else id) nnoff

genNonNegOffset :: Gen Pos.Offset
genNonNegOffset = liftM (subtract 1) genPositiveOffset

genPositiveOffset :: Gen Pos.Offset
genPositiveOffset = do scale <- chooseInteger (1, 10)
                       liftM fromIntegral $ chooseInteger (1, 2^scale)
    where chooseInteger :: (Integer, Integer) -> Gen Integer
          chooseInteger = choose

genInvertibleLoc :: Gen SpLoc.SpliceLoc
genInvertibleLoc = sized $ \sz -> do ncontigs <- choose (1, sz + 1)
                                     fwdloc <- liftM (fromJust . SpLoc.fromContigs) 
                                               $ genContigs ncontigs
                                     rc <- arbitrary
                                     if rc then return $ revCompl fwdloc else return fwdloc
    where genContigs = liftM (reverse . foldl' intervalsToContigs []) . genIntervals
          genIntervals nints = replicateM nints $ liftM2 (,) genPositiveOffset genPositiveOffset
          intervalsToContigs [] (init5, len) = [Loc.fromPosLen (Pos.Pos init5 Plus) len]
          intervalsToContigs prevs@(prev:_) (nextoffset, nextlen)
              = let !prevend = Loc.offset5 prev + Loc.length prev
                in (Loc.fromPosLen (Pos.Pos (prevend + nextoffset) Plus) nextlen):prevs

instance Arbitrary SpLoc.SpliceLoc where
  arbitrary = genInvertibleLoc

genTranscript :: Gen Transcript
genTranscript = do gene <- arbitrary
                   trx <- arbitrary
                   seqloc <- OnSeq <$> arbitrary <*> genInvertibleLoc
                   let len = fromIntegral . Loc.length . unOnSeq $ seqloc
                   s <- Pos.Offset <$> choose (0, len - 1)
                   e <- Pos.Offset <$> choose (fromIntegral s, len - 1)
                   let cds = Loc.fromStartEnd s e
                   return $! Transcript { geneId = gene, trxId = trx, location = seqloc, cds = Just $! cds }

instance Show Transcript where
  show t = "Transcript " ++ (show . BS.unpack . unSeqLabel . geneId $ t)

data Test = forall t . Testable t => T String t

runTest :: Test -> IO ()
runTest (T name test) = do
    putStr $ name ++ replicate (40 - length name) '.' ++ "  "
    quickCheckWith args test
      where args = stdArgs -- { maxDiscard = 100000 }

