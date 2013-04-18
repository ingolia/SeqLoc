{-# LANGUAGE OverloadedStrings #-}
module Main
       where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import System.Environment
import System.IO

import System.Console.CmdTheLine

import Bio.SeqLoc.GTF
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Transcript
import Bio.SeqLoc.TranscriptTable

main :: IO ()
main = run ( gtfIntrons, info )
  where info = defTI { termName = "gtf-introns"
                     , version = "0.0"
                     , termDoc = "Generate GTF file of introns from GTF file of transcripts"
                     }
        gtfIntrons = introns <$> naming <*> gtfin
        introns nam = readGtfTranscripts >=> mapM_ BS.putStr . intronGtf
          where intronGtf = transcriptsToGtf . concatMap (transcriptIntrons nam)
        
transcriptsToGtf :: [Transcript] -> [BS.ByteString]
transcriptsToGtf = map (transcriptToGtf "exons-to-introns")
        
transcriptIntrons :: Naming -> Transcript -> [Transcript]
transcriptIntrons nam trx = zipWith intronTranscript [1..] . junctions $ sploc
  where (OnSeq refname sploc) = location trx
        intronTranscript idx jct = Transcript geneid trxid loc Nothing
          where geneid = intronGeneId nam idx trx -- toSeqLabel . flip BS.append "_introns" . unSeqLabel . geneId $ trx
                trxid = intronTrxId nam idx trx -- toSeqLabel . flip BS.append trxsuffix . unSeqLabel . trxId $ trx
                loc = OnSeq refname (fromMaybe noLoc $ SpLoc.fromContigs [ intron jct ])
                noLoc = error $ "Unable to create singleton SpLoc from " ++ (BS.unpack . repr) jct
                -- trxsuffix = "_intron" `BS.append` (BS.pack . show $ idx)                
                
naming :: Term Naming
naming = value $ vFlag SameGene 
         [ ( DifferentGenes, (optInfo [ "different-genes" ]) { optDoc = "Different GTF genes for each intron" } )
         , ( SameGene, (optInfo [ "same-gene" ]) { optDoc = "Different GTF transcripts (but the same GTF gene) for each intron" })
         , ( SameTranscript, (optInfo [ "same-transcript" ]) { optDoc = "Same GTF transcript for each intron" } )
         ]

gtfin :: Term String
gtfin = required $ pos 0 Nothing $ posInfo { posName = "GTF", posDoc = "Input GTF file" }

data Naming = DifferentGenes
            | SameGene
            | SameTranscript
              
intronGeneId :: Naming -> Int -> Transcript -> SeqLabel
intronGeneId naming idx trx = toSeqLabel . flip BS.append genesuffix . unSeqLabel . geneId $ trx
  where genesuffix = case naming of
          DifferentGenes -> "_intron" `BS.append` (BS.pack . show $ idx)
          _ -> "_introns"
          
intronTrxId :: Naming -> Int -> Transcript -> SeqLabel
intronTrxId naming idx trx = toSeqLabel . flip BS.append trxsuffix . unSeqLabel . geneId $ trx
  where trxsuffix = case naming of
          SameTranscript -> "_introns"
          _ -> "_intron" `BS.append` (BS.pack . show $ idx)