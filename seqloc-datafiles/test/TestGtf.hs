{-# LANGUAGE OverloadedStrings #-}
module Main
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import System.IO

import Bio.SeqLoc.GTF
import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript
import Bio.SeqLoc.TranscriptTable

main :: IO ()
main = do trx <- readGtfTranscripts "/data/genomes/Homo_sapiens/hg19_knownGene.gtf"
          writeTable "test/out.txt" trx
          trx' <- readTable "test/out.txt"
          return ()
          