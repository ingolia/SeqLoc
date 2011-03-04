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

main :: IO ()
main = getArgs >>= handleOpt . getOpt RequireOrder optDescrs
    where handleOpt (_,    _,         errs@(_:_)) = usage (unlines errs)
          handleOpt (args, [gtf], []) = either usage (doGtfToBed gtf) $ argsToConf args
          handleOpt (_,    _,     []) = usage "Specify exactly one GTF file"
          usage errs = do prog <- getProgName
                          hPutStr stderr $ usageInfo prog optDescrs
                          hPutStrLn stderr errs

doGtfToBed :: FilePath -> Conf -> IO ()
doGtfToBed gtf conf = readGtfTranscripts gtf >>= writeTranscriptsOut
  where writeTranscriptsOut trxs = withOutHandle conf $ \h -> 
          mapM_ (BS.hPutStrLn h . transcriptToBedStd) trxs

withOutHandle :: Conf -> (Handle -> IO a) -> IO a
withOutHandle conf m = maybe (m stdout) (\outname -> withFile outname WriteMode m) $ confOutput conf

data Conf = Conf { confOutput :: !(Maybe FilePath) 
                 } deriving (Show)

data Arg = ArgOutput { unArgOutput :: !String }
         deriving (Show, Read, Eq, Ord)

argOutput :: Arg -> Maybe String
argOutput (ArgOutput del) = Just del
argOutput _ = Nothing

optDescrs :: [OptDescr Arg]
optDescrs = [ Option ['o'] ["output"] (ReqArg ArgOutput "OUTFILE") "Output filename"
            ]

argsToConf :: [Arg] -> Either String Conf
argsToConf = runReaderT conf
    where conf = Conf <$> 
                 findOutput
          findOutput = ReaderT $ return . listToMaybe . mapMaybe argOutput
