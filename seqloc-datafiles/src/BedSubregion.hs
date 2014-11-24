{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
module Main
       where  

import Control.Applicative
import Control.Exception
import qualified Control.Exception.Lifted as E
import Control.Monad
import Control.Monad.Base
import Control.Monad.IO.Class
import qualified Control.Monad.Trans.Resource as R
import qualified Data.ByteString.Char8 as BS
import qualified Data.Conduit as C
import qualified Data.Conduit.Binary as CB
import qualified Data.Conduit.List as C
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
                  where writeSubregion sl = let t' = t { location = (location t) { unOnSeq = sl }, cds = Nothing }
                                            in liftIO $ BS.hPutStrLn hout $ Bed.transcriptToBedStd t'
        
data TranscriptRegion = WholeTrx | Utr5 | Cds | Utr3 deriving (Show, Read, Eq, Ord, Bounded, Enum)

trxRegion :: TranscriptRegion -> Transcript -> Maybe Loc.ContigLoc
trxRegion WholeTrx trx = let sploc = unOnSeq . location $ trx
                         in Just $! Loc.fromPosLen (Pos.Pos 0 Plus) (Loc.length sploc)
trxRegion Utr5 trx = utr5 trx
trxRegion Utr3 trx = utr3 trx
trxRegion Cds trx = cds trx
                         
data TooLong = TooLongExtend | TooLongTruncate | TooLongDiscard deriving (Show, Read, Eq, Ord, Bounded, Enum)

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

validRegionSpec :: RegionSpec -> Bool
validRegionSpec (RegionSpec _rgn (Just startoff) (Just len) Nothing       _toolong) = True
validRegionSpec (RegionSpec _rgn (Just startoff) Nothing    (Just endoff) _toolong) = True
validRegionSpec (RegionSpec _rgn Nothing         (Just len) (Just endoff) _toolong) = True
validRegionSpec _ = False

regionSpliceLoc :: RegionSpec -> Transcript -> Maybe Loc.SpliceLoc
regionSpliceLoc (RegionSpec rgn (Just startoff) (Just len) Nothing toolong) trx
  = do base <- trxRegion rgn trx
       cloc <- handleEnds toolong base $ Loc.fromPosLen (Pos.slide (Loc.startPos base) startoff) len
       return $! clocOutofExtended cloc (unOnSeq . location $ trx)
regionSpliceLoc rs _ = error $ "Unimplemented region selection " ++ show rs

data Conf = Conf { cInput :: !FilePath
                 , cOutput :: !(Maybe FilePath)
                 , cRegionSpec :: !RegionSpec
                 }

cOutFile :: Conf -> FilePath
cOutFile conf = fromMaybe defaultOutput . cOutput $ conf
  where defaultOutput = (dropExtension . cInput $ conf) ++ "_subregion" ++ (takeExtension . cInput $ conf)

argConf :: Term Conf
argConf = Conf <$> bedin <*> bedout <*> regionspec

bedin :: Term FilePath
bedin = required $ opt Nothing $ ( optInfo [ "i" ])
  { optName = "INPUT.BED", optDoc = "BED input" }

bedout :: Term (Maybe FilePath)
bedout = value $ opt Nothing $ ( optInfo [ "o" ])
  { optName = "OUTPUT.BED", optDoc = "BED output" }

relregion :: Term TranscriptRegion
relregion = ret . fmap validate . value $
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

toolong :: Term TooLong
toolong = ret . fmap validate . value $
          vFlag Nothing $
          [ (Just TooLongExtend,   (optInfo [ "extend" ])   { optDoc = "Extend beyond reference region" })
          , (Just TooLongTruncate, (optInfo [ "truncate" ]) { optDoc = "Truncate to lie within reference region " })
          , (Just TooLongDiscard,  (optInfo [ "discard" ])  { optDoc = "Discard when lying outside reference region" })
          ]
  where validate :: Maybe TooLong -> Err TooLong
        validate = maybe noarg return
        noarg = msgFail . PP.text $
                "Specify how subregions extending outside the reference region should be handled (truncation etc.)"

startoff :: Term (Maybe Int)
startoff = value $ opt Nothing $ ( optInfo [ "start-off" ]) { optName = "DELTA-START", optDoc = "Offset of start position, positive is more 3'" }

endoff :: Term (Maybe Int)
endoff = value $ opt Nothing $ ( optInfo [ "end-off" ]) { optName = "DELTA-END", optDoc = "Offset of end position, positive is more 3'" }

rgnlen :: Term (Maybe Int)
rgnlen = value $ opt Nothing $ ( optInfo [ "length" ]) { optName = "LENGTH", optDoc = "Length of the sub-region" }

regionspec :: Term RegionSpec
regionspec = ret . fmap validate $
             RegionSpec <$>
             relregion <*>
             (liftM fromIntegral <$> startoff) <*>
             (liftM fromIntegral <$> rgnlen) <*>
             (liftM fromIntegral <$> endoff) <*>
             toolong
  where validate :: RegionSpec -> Err RegionSpec
        validate rs = if validRegionSpec rs
                      then return rs
                      else msgFail . PP.text $
                           "Specify exactly two of start offset, length, and end offset"

throwerr :: (MonadBase IO m) => String -> m a
throwerr = E.ioError . userError
