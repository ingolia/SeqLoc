module Main
       where

import Control.Applicative
import Control.Monad.Reader
import qualified Data.ByteString.Char8 as BS
import Data.Maybe

import System.Console.GetOpt
import System.Environment
import System.IO

import Bio.SeqLoc.Bed
import Bio.SeqLoc.GTF
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.LocRepr
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [gtf], []) = either usage (doGtfToBed gtf) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify exactly one GTF file"
          usage errs = do prog <- getProgName
                          hPutStr stderr $ usageInfo prog optDescrs
                          hPutStrLn stderr errs

doGtfToBed :: FilePath -> Conf -> IO ()
doGtfToBed gtf conf = readGtfTranscriptsErr gtf >>= handleErrors >>= writeTranscriptsOut . map adjustStop
  where writeTranscriptsOut trxs = withOutHandle conf $ \h -> 
          mapM_ (BS.hPutStrLn h . transcriptToBedStd) trxs
        adjustStop = if confStopIncluded conf then removeDoubleStop else id
        removeDoubleStop trx = trx { cds = cds trx >>= shortenLoc }
        shortenLoc l = Loc.clocOutof (Loc.fromBoundsStrand 0 end Plus) l
          where end = max 0 $ Loc.length l - 4
        handleErrors (trxs, errs) = do maybe putErrs filePutErrs $ confBadTranscripts conf
                                       return trxs
          where putErrs = hPutErrs stderr
                filePutErrs badtx = withFile badtx WriteMode hPutErrs
                hPutErrs h = hPutStr h (unlines errs)

withOutHandle :: Conf -> (Handle -> IO a) -> IO a
withOutHandle conf m = maybe (m stdout) (\outname -> withFile outname WriteMode m) $ confOutput conf

data Conf = Conf { confOutput :: !(Maybe FilePath) 
                 , confBadTranscripts :: !(Maybe FilePath)
                 , confStopIncluded :: !Bool
                 } deriving (Show)

data Arg = ArgOutput { unArgOutput :: !String }
         | ArgBadTranscripts { unArgBadTranscripts :: !String }
         | ArgStopIncluded
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput del) = Just del
argOutput _ = Nothing

argBadTranscripts :: Arg -> Maybe String
argBadTranscripts (ArgBadTranscripts bad) = Just bad
argBadTranscripts _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"]        (ReqArg ArgOutput "OUTFILE") "Output filename"
            , Option ['b'] ["bad-transcripts"] (ReqArg ArgBadTranscripts "BADFILE") "Write bad transcript list to file"
            , Option ['s'] ["stop-included"] (NoArg ArgStopIncluded)      "GTF file includes stop codon in CDS"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput <*>
                 findBadTranscripts <*>
                 (ReaderT $ return . elem ArgStopIncluded)
          findOutput = ReaderT $ return . listToMaybe . mapMaybe argOutput
          findBadTranscripts = ReaderT $ return . listToMaybe . mapMaybe argBadTranscripts