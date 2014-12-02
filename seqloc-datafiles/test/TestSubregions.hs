{-# LANGUAGE BangPatterns, ExistentialQuantification #-}
module Main
       where

import Control.Applicative
import Control.Monad
import Control.Monad.IO.Class
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import System.Directory
import System.IO
import System.Process
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

bedSubregion = "./dist/build/bed-subregion/bed-subregion"

main :: IO ()
main = mapM_ runTest tests

tests :: [Test]
tests = [ T "Write random bed"          test_writeRandomBed
        , T "Whole transcript"          test_wholeTranscript
        ]

test_writeRandomBed :: Property
test_writeRandomBed = monadicIO $
                      forAllM (listOf1 genTranscript) $ \trxs ->
                      run . BS.writeFile "tmpf.bed" . BS.unlines . map Bed.transcriptToBedStd $ trxs

test_wholeTranscript :: Property
test_wholeTranscript = test_subtranscript (Just . stripCDS) [ "--whole-trx", "--discard", "--start=0", "--end=0" ]
  where stripCDS t = t { cds = Nothing }

test_subtranscript :: (Transcript -> Maybe Transcript) -> [String] -> Property
test_subtranscript subtrx args = monadicIO $
                                 forAllM (listOf1 genTranscript) $ \trxs ->
                                 let subtrxs = mapMaybe subtrx trxs
                                 in run $ do (fullname, hfull) <- openTempFile "." "test-subtrx-full.bed"
                                             (subname, hsub) <- openTempFile "." "test-subtrx-sub.bed"
                                             (outname, hout) <- openTempFile "." "test-subtrx-out.bed"
                                             hClose hout
                                             removeFile outname
                                             BS.hPutStr hfull . BS.unlines . map Bed.transcriptToBedStd $ trxs
                                             hClose hfull
                                             BS.hPutStr hsub . BS.unlines . map Bed.transcriptToBedStd $ subtrxs
                                             hClose hsub
                                             callProcess bedSubregion (args ++ [ "-i", fullname, "-o", outname ])
                                             callProcess "diff" [ subname, outname ]
                                             removeFile fullname
                                             removeFile subname
                                             removeFile outname

--test_cdsStartLength :: Property
--test_cdsStartLength = monadicIO $
--                      forAllM (listOf1 genTranscript) $ \trxs ->
--                      let subtrxs = su

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
                   cds <- arbitrary >>= \hascds ->
                     if hascds
                        then let len = fromIntegral . Loc.length . unOnSeq $ seqloc
                             in do s <- Pos.Offset <$> choose (0, len - 1)
                                   e <- Pos.Offset <$> choose (fromIntegral s, len - 1)
                                   return . Just $! Loc.fromStartEnd s e
                        else return Nothing
                   return $! Transcript { geneId = gene, trxId = trx, location = seqloc, cds = cds }

instance Show Transcript where
  show t = "Transcript " ++ (show . BS.unpack . unSeqLabel . geneId $ t)

data Test = forall t . Testable t => T String t

runTest :: Test -> IO ()
runTest (T name test) = do
    putStr $ name ++ replicate (40 - length name) '.' ++ "  "
    quickCheckWith args test
      where args = stdArgs -- { maxDiscard = 100000 }

