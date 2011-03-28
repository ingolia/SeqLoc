{-| Minimal tab-delimited annotation of 'Transcript' locations. Each
'Transcript' has one line, containing the 'geneId' and 'trxId' fields,
followed by the 'LocRepr' representation of the 'SpliceSeqLoc'
location of the transcript, and then the location of the CDS within
the transcript or "n/a" if there is no CDS. -}

module Bio.SeqLoc.TranscriptTable
       ( readTable, parseLine, writeTable, unparseLine )
       where 

import Control.Monad
import qualified Data.ByteString.Char8 as BS

import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.OnSeq
import Bio.SeqLoc.Transcript

-- | Read a transcript table file as a list of annotated transcripts
readTable :: FilePath -> IO [Transcript]
readTable = BS.readFile >=> mapM parseLineErr . BS.lines
  where parseLineErr l = maybe err return $! parseLine l
          where err = ioError . userError $ "Bad transcript table line " ++ show l

-- | Parse a transcript table line, not including the newline
parseLine :: BS.ByteString -> Maybe Transcript
parseLine l = case BS.split '\t' l of
  [ geneid, trxid, trxlocstr, cdslocstr ] -> do
    trxloc <- unreprMaybe trxlocstr
    cdsloc <- if cdslocstr == na
                 then return Nothing
                 else liftM Just $! unreprMaybe cdslocstr
    return $! Transcript (SeqName . BS.copy $ geneid) (SeqName . BS.copy $ trxid) trxloc cdsloc
  _ -> Nothing

-- | Write a transcript table file
writeTable :: FilePath -> [Transcript] -> IO ()
writeTable outname = BS.writeFile outname . BS.unlines . map unparseLine

-- | Produce a single transcript table line for a 'Transcript', newline not included.
unparseLine :: Transcript -> BS.ByteString
unparseLine trx = BS.intercalate (BS.singleton '\t') $ [ unSeqName . geneId $ trx
                                                       , unSeqName . trxId $ trx
                                                       , repr . location $ trx
                                                       , maybe na repr . cds $ trx
                                                       ]
                  
na :: BS.ByteString
na = BS.pack "n/a"