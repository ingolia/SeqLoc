{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
  where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.Char
import qualified Data.HashMap.Strict as HM
import Data.Maybe
import Numeric
import System.Exit
import System.FilePath
import System.IO

import qualified Data.Attoparsec.Zepto as ZP
import qualified Data.Vector.Unboxed as U
import System.Console.CmdTheLine

import qualified Control.Monad.Trans.Resource as R
import qualified Data.Conduit as C
import qualified Data.Conduit.Binary as CB
import qualified Data.Conduit.List as C

import qualified Bio.SeqLoc.Bed as Bed
import qualified Bio.SeqLoc.LocMap as LM
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript
import Bio.SeqLoc.ZeptoUtils

main :: IO ()
main = run ( g2t, info )
  where info = defTI { termName = "genome-to-trx"
                     , version = "0.0"
                     , termDoc = "Convert genome coordinates to transcriptome coordinates"
                     }
        g2t = genomeToTrx <$> trxBed <*> coordInput <*> outputFile <*> argConf

genomeToTrx :: String -> String -> Maybe String -> Conf -> IO ()
genomeToTrx trx input moutput conf
  = do trxs <- liftM (LM.transcriptSeqLocMap 100000) $ Bed.readBedTranscripts trx
       R.runResourceT $ CB.sourceFile input C.$= remapCoordLines trxs conf C.$$ CB.sinkFile output
  where output = fromMaybe defaultOutput moutput
        defaultOutput = let (base, ext) = splitExtension input
                        in (base ++ "_trx") <.> ext

data InputType = InputBed | InputBedPE deriving (Read, Show, Ord, Eq)
data Strandedness = FwdOnly | RevOnly | Both deriving (Read, Show, Ord, Eq)

data Conf = Conf { cInputType :: !InputType,
                   cStrandedness :: !Strandedness,
                   cReportNoHit :: !Bool
                 } deriving (Read, Show)

trxBed :: Term String
trxBed = required $ opt Nothing $ (optInfo [ "t", "transcripts" ])
  { optName = "TRANSCRIPTS.BED", optDoc = "Transcript BED file" }

coordInput :: Term String
coordInput = required $ opt Nothing $ (optInfo [ "i", "input" ])
  { optName = "INPUT.TXT", optDoc = "Input coordinates" }

outputFile :: Term (Maybe String)
outputFile = value $ opt Nothing $ (optInfo [ "o", "output" ])
  { optName = "OUTPUT.TXT", optDoc = "Output coordinates" }

argConf :: Term Conf
argConf = Conf <$>
          confInputType <*>
          confStrandedness <*>
          confReportNoHit
  where confInputType :: Term InputType
        confInputType = value $ vFlag InputBed [(InputBedPE, (optInfo [ "bedpe" ]) { optName = "Bed PE (bedtools) format", optDoc = "Bed PE (bedtools) format" })]
        confStrandedness :: Term Strandedness
        confStrandedness = value $ vFlag FwdOnly
                           [( FwdOnly, (optInfo [ "fwd"  ]) { optName = "Forward feature strand only", optDoc = "Forward feature strand only" }),
                            ( RevOnly, (optInfo [ "rev"  ]) { optName = "Reverse feature strand only", optDoc = "Reverse feature strand only" }),
                            ( Both,    (optInfo [ "both" ]) { optName = "Both feature strands",        optDoc = "Both feature strands" })]
        confReportNoHit :: Term Bool
        confReportNoHit = value $ flag $ (optInfo [ "n", "no-hit" ]) { optName = "Report no-hit lines", optDoc = "Report no-hit lines" }

remapCoordLines :: (Monad m) => LM.SeqLocMap Transcript -> Conf -> C.Conduit BS.ByteString m BS.ByteString
remapCoordLines trxs conf = CB.lines C.=$= C.concatMapM remapCoordLine C.=$= C.map (flip BS.append "\n")
  where remapCoordLine :: (Monad m) => BS.ByteString -> m [BS.ByteString]
        remapCoordLine = case cInputType conf of
          InputBed -> remapBedLine trxs conf
          etc -> fail $ "Input type " ++ show etc ++ " not implemented"

remapBedLine :: (Monad m) => LM.SeqLocMap Transcript -> Conf -> BS.ByteString -> m [BS.ByteString]
remapBedLine trxs conf l = case BS.split '\t' l of
  (chr:startBS:endBS:name:score:strandBS:rest)
    -> do start <- either (const . fail $ "Bad start in " ++ show l) return $ ZP.parse decimal startBS
          end <- either (const . fail $ "Bad end in " ++ show l) return $ ZP.parse decimal endBS
          strand <- case strandBS of
            "+" -> return Plus
            "-" -> return Minus
            _ -> fail $ "Bad strand " ++ show strandBS ++ " in " ++ show l
          let gloc = OnSeq (toSeqLabel chr) (Loc.fromBoundsStrand start (end - 1) strand)
          let tlocs = remapLoc trxs conf gloc
              bedLine tsloc = BS.intercalate "\t" [ unSeqLabel . onSeqLabel $ tsloc
                                                  , BS.pack . show . Pos.unOff . fst . Loc.bounds . unOnSeq $ tsloc
                                                  , BS.pack . show . (+ 1) . Pos.unOff . snd . Loc.bounds . unOnSeq $ tsloc
                                                  , name, score
                                                  , if ((Loc.strand . unOnSeq) tsloc == Plus) then "+" else "-" ]
          if null tlocs && cReportNoHit conf
             then return $! [BS.intercalate "\t" $ [ "N/A", ".", ".", name, score, "." ]]
             else return $! map bedLine tlocs

remapLoc :: LM.SeqLocMap Transcript -> Conf -> ContigSeqLoc -> [ContigSeqLoc]
remapLoc trxs conf l = mapMaybe locInto cands
  where candList = LM.querySeqLoc l trxs
        cands = HM.elems . HM.fromList . map (\t -> (trxId t, t)) $ candList
        locInto t = let tsloc = location t
                    in Loc.clocInto (unOnSeq l) (unOnSeq tsloc) >>= \tloc ->
                    if (Loc.strand tloc == Plus && cStrandedness conf == RevOnly) ||
                       (Loc.strand tloc == Minus && cStrandedness conf == FwdOnly)
                       then Nothing
                       else Just (OnSeq (trxId t) tloc)
                   
                          
