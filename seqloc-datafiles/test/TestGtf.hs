{-# LANGUAGE OverloadedStrings #-}
module Main
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import System.IO

import Bio.SeqLoc.GTF
import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript
import Bio.SeqLoc.TranscriptTable

main :: IO ()
main = do trx <- readGtfTranscripts "/data/genomes/Homo_sapiens/hg19_knownGene.gtf"
          writeTable "test/gtf-out.txt" trx
          trx' <- readTable "test/gtf-out.txt"
          let trx10 = take 10 trx'
          BS.writeFile "test/gtf-out10.gtf" . BS.concat . map (transcriptToGtf "TestGtf") $ trx10
          trx10' <- readGtfTranscripts "test/gtf-out10.gtf"
          print $ (sort . map location $ trx10) == (sort . map location $ trx10')
          return ()
          