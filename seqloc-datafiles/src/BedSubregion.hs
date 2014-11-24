{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
module Main
       where  

import Control.Applicative
import qualified Control.Exception.Lifted as E
import Control.Monad
import Control.Monad.Base
import Control.Monad.IO.Class
import qualified Control.Monad.Trans.Resource as R
import qualified Data.ByteString.Char8 as BS
import qualified Data.Conduit as C
import qualified Data.Conduit.Binary as CB
import qualified Data.Conduit.List as C
import Data.List
import Data.Maybe

import qualified Bio.SeqLoc.Bed as Bed
import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as Loc
import Bio.SeqLoc.Strand
import Bio.SeqLoc.Transcript

import System.Console.CmdTheLine
import System.FilePath
import System.IO
import qualified Text.PrettyPrint as PP

main :: IO ()
main = run ( bedsub, info )
  where info = defTI { termName = "bed-subregion"
                     , version = "2014-11-23"
                     , termDoc = "Process BED file to extract specific sub-regions"
                     }
        bedsub = bedSubregions <$> argConf

bedSubregions :: Conf -> IO ()
bedSubregions conf = R.runResourceT $
                     CB.sourceFile (cInput conf) C.$=
                     Bed.bedConduit C.$$
                     handleTranscripts conf

handleTranscripts :: (MonadIO m, R.MonadResource m) => Conf -> C.Sink Transcript m ()
handleTranscripts conf = C.bracketP (openFile (cOutFile conf) WriteMode) hClose loop
  where loop hout = C.head >>= maybe (return ()) (\t -> handleOne t >> loop hout)
          where handleOne t = maybe (return ()) writeSubregion $ regionSpliceLoc (cRegionSpec conf) t
                  where writeSubregion sl = liftIO . BS.hPutStrLn hout . Bed.transcriptToBedStd $ subtranscript sl
                        subtranscript sl = Transcript { geneId = gene', trxId = trx', location = loc', cds = Nothing }
                          where loc' = (location t) { unOnSeq = sl }
                                gene' = toSeqLabel . flip BS.append suffix . unSeqLabel . geneId $ t
                                trx' = toSeqLabel . flip BS.append suffix . unSeqLabel . trxId $ t
                                suffix = if cRenameFeatures conf
                                            then BS.empty
                                            else BS.pack . ('_' :) . regionSpecName . cRegionSpec $ conf
        
data TranscriptRegion = WholeTrx | Utr5 | Cds | Utr3 deriving (Show, Read, Eq, Ord, Bounded, Enum)

regionName :: TranscriptRegion -> String
regionName WholeTrx = "trx"
regionName Utr5 = "utr5"
regionName Cds = "cds"
regionName Utr3 = "utr3"

trxRegion :: TranscriptRegion -> Transcript -> Maybe Loc.ContigLoc
trxRegion WholeTrx trx = let sploc = unOnSeq . location $ trx
                         in Just $! Loc.fromPosLen (Pos.Pos 0 Plus) (Loc.length sploc)
trxRegion Utr5 trx = utr5 trx
trxRegion Utr3 trx = utr3 trx
trxRegion Cds trx = cds trx
                         
data TooLong = TooLongExtend | TooLongTruncate | TooLongDiscard deriving (Show, Read, Eq, Ord, Bounded, Enum)

tooLongName :: TooLong -> String
tooLongName TooLongExtend = "ext"
tooLongName TooLongTruncate = "trunc"
tooLongName TooLongDiscard = "disc"

handleEnds :: TooLong -> Loc.ContigLoc -> Loc.ContigLoc -> Maybe Loc.ContigLoc
handleEnds TooLongExtend _base cloc = Just cloc
handleEnds TooLongTruncate base cloc = let start = max (Pos.offset . Loc.startPos $ base) (Pos.offset . Loc.startPos $ cloc)
                                           end = min (Pos.offset . Loc.endPos $ base) (Pos.offset . Loc.endPos $ cloc)
                                       in if start < end
                                          then Just $! Loc.fromStartEnd start end
                                          else Nothing
handleEnds TooLongDiscard base cloc = if ((Pos.offset . Loc.startPos $ base) <= (Pos.offset . Loc.startPos $ cloc)
                                          && (Pos.offset . Loc.endPos $ base) >= (Pos.offset . Loc.endPos $ cloc))
                                      then Just cloc
                                      else Nothing

clocOutofExtended :: (Loc.Location l, LocRepr l) => Loc.ContigLoc -> l -> l
clocOutofExtended cloc l0 = fromMaybe badExtension $ Loc.clocOutof (Loc.slide startxtn cloc) lext
  where startxtn = negate . Pos.offset . Loc.startPos $ cloc
        endxtn = (Pos.offset . Loc.endPos $ cloc) + 1 - (Loc.length l0 - 1)
        lext = Loc.extend (startxtn, endxtn) l0
        badExtension = error . BS.unpack $ BS.unwords [ "Bad extension: ", repr cloc, " in ", repr l0, " to ", repr lext ]

data RegionSpec = RegionSpec { rsRegion :: !TranscriptRegion
                             , rsStartOffset :: !(Maybe Pos.Offset)
                             , rsLength :: !(Maybe Pos.Offset)
                             , rsEndOffset :: !(Maybe Pos.Offset)
                             , rsTooLong :: !TooLong
                             }
                deriving (Show)

regionSpecName :: RegionSpec -> String
regionSpecName (RegionSpec rgn (Just startoff) (Just len) Nothing toolong)
  = intercalate "_" [ regionName rgn, "start" ++ showSigned startoff, "length" ++ (show . Pos.unOff $ len), tooLongName toolong ]
regionSpecName (RegionSpec rgn (Just startoff) Nothing (Just endoff) toolong)
  = intercalate "_" [ regionName rgn, "start" ++ showSigned startoff, "end" ++ show endoff, tooLongName toolong ]
regionSpecName (RegionSpec rgn Nothing (Just len) (Just endoff) toolong)
  = intercalate "_" [ regionName rgn, "length" ++ (show . Pos.unOff $ len), "end" ++ show endoff, tooLongName toolong ]
regionSpecName rs = error $ "regionSpecName: invalid region specification " ++ show rs

showSigned :: Pos.Offset -> String
showSigned z = case show . Pos.unOff $ z of
  [] -> []
  str@('-':_) -> str
  str -> '+' : str

validRegionSpec :: RegionSpec -> Bool
validRegionSpec (RegionSpec _rgn (Just _startoff) (Just _len) Nothing        _toolong) = True
validRegionSpec (RegionSpec _rgn (Just _startoff) Nothing     (Just _endoff) _toolong) = True
validRegionSpec (RegionSpec _rgn Nothing          (Just _len) (Just _endoff) _toolong) = True
validRegionSpec _ = False

regionSpliceLoc :: RegionSpec -> Transcript -> Maybe Loc.SpliceLoc
regionSpliceLoc (RegionSpec rgn (Just startoff) (Just len) Nothing toolong) trx
  = do base <- trxRegion rgn trx
       cloc <- handleEnds toolong base $ Loc.fromPosLen (Pos.slide (Loc.startPos base) startoff) len
       return $! clocOutofExtended cloc (unOnSeq . location $ trx)
regionSpliceLoc (RegionSpec rgn (Just startoff) Nothing (Just endoff) toolong) trx
  = do base <- trxRegion rgn trx
       let start = (Pos.offset . Loc.startPos $ base) + startoff
           end = (Pos.offset . Loc.endPos $ base) + endoff
       if start <= end
          then do cloc <- handleEnds toolong base $ Loc.fromStartEnd start end
                  return $! clocOutofExtended cloc (unOnSeq . location $ trx)
          else Nothing
regionSpliceLoc (RegionSpec rgn Nothing (Just len) (Just endoff) toolong) trx
  = do base <- trxRegion rgn trx
       let end = (Pos.offset . Loc.endPos $ base) + endoff
           start = 1 + end - len
       if start <= end
          then do cloc <- handleEnds toolong base $ Loc.fromStartEnd start end
                  return $! clocOutofExtended cloc (unOnSeq . location $ trx)
          else Nothing
regionSpliceLoc rs _ = error $ "Invalid region selection " ++ show rs

data Conf = Conf { cInput :: !FilePath
                 , cOutput :: !(Maybe FilePath)
                 , cRegionSpec :: !RegionSpec
                 , cRenameFeatures :: !Bool
                 }

cOutFile :: Conf -> FilePath
cOutFile conf = fromMaybe defaultOutput . cOutput $ conf
  where defaultOutput = (dropExtension . cInput $ conf) ++ "_" ++ regionSpecName (cRegionSpec conf) ++ (takeExtension . cInput $ conf)

argConf :: Term Conf
argConf = Conf <$> argInput <*> argOutput <*> regionspec <*> argRename

argInput :: Term FilePath
argInput = required $ opt Nothing $ ( optInfo [ "i" ])
  { optName = "INPUT.BED", optDoc = "BED input" }

argOutput :: Term (Maybe FilePath)
argOutput = value $ opt Nothing $ ( optInfo [ "o" ])
  { optName = "OUTPUT.BED", optDoc = "BED output" }

argRegion :: Term TranscriptRegion
argRegion = ret . fmap validate . value $
            vFlag Nothing $
            [ (Just WholeTrx, (optInfo [ "whole-trx" ]) { optDoc = "Region relative to whole transcript" })
            , (Just Utr5,     (optInfo [ "utr5" ])      { optDoc = "Region relative to 5' UTR" })
            , (Just Utr3,     (optInfo [ "utr3" ])      { optDoc = "Region relative to 3' UTR" })
            , (Just Cds,      (optInfo [ "cds" ])       { optDoc = "Region relative to CDS" })
            ]
  where validate :: Maybe TranscriptRegion -> Err TranscriptRegion
        validate = maybe noRegion return
        noRegion = msgFail . PP.text $
                   "Specify a reference transcription region (whole transcript, CDS, etc.)"

argTooLong :: Term TooLong
argTooLong = ret . fmap validate . value $
             vFlag Nothing $
             [ (Just TooLongExtend,   (optInfo [ "extend" ])   { optDoc = "Extend beyond reference region" })
             , (Just TooLongTruncate, (optInfo [ "truncate" ]) { optDoc = "Truncate to lie within reference region " })
             , (Just TooLongDiscard,  (optInfo [ "discard" ])  { optDoc = "Discard when lying outside reference region" })
             ]
  where validate :: Maybe TooLong -> Err TooLong
        validate = maybe noarg return
        noarg = msgFail . PP.text $
                "Specify how subregions extending outside the reference region should be handled (truncation etc.)"

argStartOffset :: Term (Maybe Int)
argStartOffset = value $ opt Nothing $ ( optInfo [ "start-off" ]) { optName = "DELTA-START", optDoc = "Offset of start position, positive is more 3'" }

argEndOffset :: Term (Maybe Int)
argEndOffset = value $ opt Nothing $ ( optInfo [ "end-off" ]) { optName = "DELTA-END", optDoc = "Offset of end position, positive is more 3'" }

argLength :: Term (Maybe Int)
argLength = value $ opt Nothing $ ( optInfo [ "length" ]) { optName = "LENGTH", optDoc = "Length of the sub-region" }

regionspec :: Term RegionSpec
regionspec = ret . fmap validate $
             RegionSpec <$>
             argRegion <*>
             (liftM fromIntegral <$> argStartOffset) <*>
             (liftM fromIntegral <$> argLength) <*>
             (liftM fromIntegral <$> argEndOffset) <*>
             argTooLong
  where validate :: RegionSpec -> Err RegionSpec
        validate rs = if validRegionSpec rs
                      then return rs
                      else msgFail . PP.text $
                           "Specify exactly two of start offset, length, and end offset"

argRename :: Term Bool
argRename = value $ opt False $ ( optInfo [ "r", "rename-features" ])
            { optDoc = "Rename features to describe subregion" }

throwerr :: (MonadBase IO m) => String -> m a
throwerr = E.ioError . userError
