{-# LANGUAGE OverloadedStrings #-}
module Main
       where

import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import System.Environment
import System.IO

import Bio.SeqLoc.GTF
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Transcript
import Bio.SeqLoc.TranscriptTable

main :: IO ()
main = getArgs >>= mainWithArgs
  where mainWithArgs [gtf] = gtfIntrons gtf
        mainWithArgs _ = do prog <- getProgName
                            hPutStrLn stderr . unwords $ [ "USAGE:", prog, "<GTF>" ]
                            
gtfIntrons :: FilePath -> IO ()
gtfIntrons = readGtfTranscripts >=> mapM_ BS.putStr . intronGtf
  where intronGtf = transcriptsToGtf . concatMap transcriptIntrons
        
transcriptsToGtf :: [Transcript] -> [BS.ByteString]
transcriptsToGtf = map (transcriptToGtf "exons-to-introns")
        
transcriptIntrons :: Transcript -> [Transcript]
transcriptIntrons trx = zipWith intronTranscript [1..] . junctions $ sploc
  where (OnSeq refname sploc) = location trx
        intronTranscript idx jct = Transcript geneid trxid loc Nothing
          where geneid = toSeqLabel . flip BS.append "_introns" . unSeqLabel . geneId $ trx
                trxid = toSeqLabel . flip BS.append trxsuffix . unSeqLabel . trxId $ trx
                trxsuffix = "_intron" `BS.append` (BS.pack . show $ idx)
                loc = OnSeq refname (fromMaybe noLoc $ SpLoc.fromContigs [ intron jct ])
                noLoc = error $ "Unable to create singleton SpLoc from " ++ (BS.unpack . repr) jct