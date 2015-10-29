{-# LANGUAGE OverloadedStrings #-}
module Main
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import System.IO

import qualified Data.Iteratee as Iter

import Bio.SeqLoc.Bed
import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript
import Bio.SeqLoc.TranscriptTable

main :: IO ()
main = do trxs <- readBedTranscripts "/mnt/ingolialab/ingolia/Genomes/Saccharomyces_cerevisiae/YeastGenome/sac_cer_yassour.bed"
          hPutStrLn stderr $! "Got " ++ (show . length $ trxs) ++ " transcripts"
          trxs <- readBedTranscripts "/mnt/ingolialab/ingolia/Genomes/Homo_sapiens/GRCh38/gencode.bed"
          hPutStrLn stderr $! "Got " ++ (show . length $ trxs) ++ " transcripts"
          
-- main = do withFile "test/bed-out.txt" WriteMode $ \hout ->
--             let bedIter = bedTranscriptEnum $ Iter.mapM_ (BS.hPutStrLn hout . unparseLine)
--             in Iter.fileDriver bedIter "/data/genomes/Homo_sapiens/hg19_knownGene.bed"
--           withFile "test/bed-copy.bed" WriteMode $ \hout ->
--             let bedIter = bedTranscriptEnum $ Iter.mapM_ (BS.hPutStrLn hout . transcriptToBedStd)
--             in Iter.fileDriver bedIter "/data/genomes/Homo_sapiens/hg19_knownGene.bed"
            
          -- trx' <- readTable "test/out.txt"
          -- let trx10 = take 10 trx'
          -- BS.writeFile "test/out10.gtf" . BS.concat . map (transcriptToGtf "TestGtf") $ trx10
          -- trx10' <- readGtfTranscripts "test/out10.gtf"
          -- print $ (sort . map location $ trx10) == (sort . map location $ trx10')
          -- return ()
          
