{-# LANGUAGE ExistentialQuantification, BangPatterns, ScopedTypeVariables #-}
module Main
    where

import Control.Applicative
import Control.Monad
import qualified Data.ByteString as BSW
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as LBS
import Data.ByteString.Internal (c2w, w2c)
import Data.Char
import Data.Either
import Data.Ix (inRange)
import Data.List
import Data.Maybe
import System.Random

import Test.QuickCheck

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand
import qualified Bio.SeqLoc.SeqLike as SeqLike

import qualified Bio.SeqLoc.LocMap as LM
import qualified Bio.SeqLoc.ShiftedVector as ShV

main :: IO ()
main = mapM_ runTest tests

tests :: [Test]
tests = [ T "Strand revCompl"               test_Strand_revCompl 
        , T "Char revCompl"                 property_Char_revCompl
        , T "ByteString revCompl"           property_ByteString_revCompl
        , T "Sequence revCompl"             property_Sequence_revCompl

        , T "ShVector singleton"            property_ShVector_singleton
        , T "ShVector update1"              property_ShVector_update1
        , T "ShVector update2"              property_ShVector_update2

        , T "LocMap hit inside only"        property_LocMap_hitIn
        , T "LocMap hit all"                property_LocMap_hitAll
        , T "LocMap hit multi"              property_LocMap_hitMulti

        , T "Locatable hit within"          property_Locatable_hitWithin
          
        , T "Pos revCompl"                  test_Pos_revCompl
        , T "Pos atPos"                     property_Pos_atPos
        , T "Pos atPos2"                    property_Pos_atPos2
        , T "Pos repr"                      test_Pos_repr
          
        , T "Contig revCompl"               test_Contig_Minus
        , T "Contig pos into/outof inverse" property_ContigIntoOutof
        , T "Contig pos outof/into inverse" property_ContigOutofInto
        , T "Contig loc into/outof inverse" property_ContigLocIntoOutof
        , T "Contig loc outof/into inverse" property_ContigLocOutofInto
        , T "Contig into based on bounds"   test_Contig_IntoBounds
        , T "Contig outof based on bounds"  test_Contig_OutofBounds
        , T "Contig allPos/outof equiv"     property_Contig_allPos_outof
        , T "Contig seqData"                property_Contig_seqData
        , T "Contig seqDataPadded"          property_Contig_seqDataPadded
        , T "Contig seqData2"               property_Contig_seqData2
        , T "Contig extend/revCompl"        property_Contig_extendMinus
        , T "Contig fromStartEnd"           property_Contig_fromStartEnd
        , T "Contig fromBoundsStrand"       property_Contig_fromBoundsStrand
        , T "Contig overlaps"               property_Contig_overlaps  
        , T "Contig repr"                   test_Contig_repr
          
        , T "Loc revCompl"                  test_Loc_Minus
        , T "Loc pos into/outof inverse"    property_LocIntoOutof
        , T "Loc pos outof/into inverse"    property_LocOutofInto
        , T "Loc outof based on bounds"     test_Loc_OutofBounds
        , T "Loc loc outof/into inverse"    property_LocCLocOutofInto
        , T "Loc outof association"         property_LocOutofAssoc
        , T "Loc allPos/outof equiv"        property_Loc_allPos_outof        
        , T "Loc locOutof"                  property_SpLocOutof
        , T "Loc locOutof valid"            property_SpLocOutofGood
        , T "Loc within"                    property_Loc_Within
        , T "Loc seqData"                   property_Loc_seqData
        , T "Loc seqDataPadded"             property_Loc_seqDataPadded
        , T "SpLoc seqData2"                property_SpLoc_seqData2
        , T "SpLoc repr"                    test_SpLoc_repr
        , T "SpLoc termini/revCompl"        property_SpLoc_terminiMinus
        , T "SpLoc termini/extend"          property_SpLoc_terminiExtend

        ]


-- Bio.BioSeq.Stranded

genNtByteString :: Int -> Gen BS.ByteString
genNtByteString = liftM BS.pack . flip replicateM (elements "ACGT")

genName :: Gen SeqLabel
genName = liftM (SeqLabel . LBS.pack) $ genNameLength >>= flip replicateM genNameChar
    where genNameLength = choose (1, 20)
          genNameChar = elements $ ['a'..'z'] ++ ['A'..'Z'] ++ ['0'..'9'] ++ "-_"

instance Arbitrary SeqLabel where
    arbitrary = genName

test_revCompl :: (Eq s, Stranded s) => s -> Bool
test_revCompl s = (revCompl . revCompl) s == s

test_repr :: (LocRepr l, Eq l) => l -> Bool
test_repr l = (unreprMaybe . repr $ l) == Just l

test_Strand_revCompl :: Strand -> Bool
test_Strand_revCompl = test_revCompl

property_Char_revCompl :: Property
property_Char_revCompl = forAll (elements "ACGTacgtnN") test_revCompl

property_ByteString_revCompl :: Property
property_ByteString_revCompl = forAll (sized genNtByteString) test_revCompl

property_Sequence_revCompl :: Property
property_Sequence_revCompl
    = forAll arbitrary $ \name ->
      let mkSeq = OnSeq name
      in forAll (sized genNtByteString) $ \sequ ->
          (unOnSeq . revCompl . mkSeq) sequ == revCompl sequ

-- Bio.BioSeq.Position

test_Pos_revCompl :: Pos.Pos -> Bool
test_Pos_revCompl = test_revCompl

property_Pos_atPos :: Pos.Pos -> Property
property_Pos_atPos pos
    = forAll genPositiveOffset $ \seqlen ->
    forAll (genNtByteString $ fromIntegral seqlen) $ \sequ ->
    let actual = Pos.atPos sequ pos
    in if and [ Pos.offset pos >= 0, Pos.offset pos < seqlen ]
       then let fwdNt = BS.index sequ (fromIntegral . Pos.offset $ pos)
            in case Pos.strand pos of
              Plus ->      actual == Just fwdNt
              Minus -> actual == Just (compl $ fwdNt)
       else actual == Nothing

property_Pos_atPos2 :: Pos.Pos -> Property
property_Pos_atPos2 pos
  = forAll genPositiveOffset $ \seqlen ->
  forAll (genNtByteString $ fromIntegral seqlen) $ \sequ -> 
  and [ Pos.atPos sequ pos == Pos.atPos (LBS.fromChunks [sequ]) pos
      , Pos.atPos sequ pos == Pos.atPos (BS.unpack sequ) pos
      ]

test_Pos_repr :: Pos.Pos -> Bool
test_Pos_repr = test_repr

-- Bio.BioSeq.Location

instance Arbitrary Strand where
    arbitrary = elements [Plus, Minus]

instance Arbitrary Pos.Pos where
    arbitrary = liftM2 Pos.Pos genOffset arbitrary

instance Arbitrary Loc.ContigLoc where
    arbitrary = liftM2 Loc.fromPosLen arbitrary genPositiveOffset

test_Contig_Minus :: Loc.ContigLoc -> Bool
test_Contig_Minus = test_revCompl

property_ContigIntoOutof :: Loc.ContigLoc -> Pos.Pos -> Property
property_ContigIntoOutof contig pos
    = let !mInpos = Loc.posInto pos contig
          !mOutpos = mInpos >>= flip Loc.posOutof contig                     
      in (isJust mInpos) ==> mOutpos == (Just pos)

property_ContigOutofInto :: Loc.ContigLoc -> Pos.Pos -> Property
property_ContigOutofInto contig pos
    = let !mOutpos = Loc.posOutof pos contig
          !mInpos = mOutpos >>= flip Loc.posInto contig
      in (isJust mOutpos) ==> mInpos == (Just pos)

property_ContigLocIntoOutof :: Loc.ContigLoc -> Loc.ContigLoc -> Property
property_ContigLocIntoOutof subcloc supercloc
    = let !mIncloc = Loc.clocInto subcloc supercloc
          !mOutcloc = mIncloc >>= flip Loc.clocOutof supercloc
      in (isJust mIncloc) ==> mOutcloc == (Just subcloc)

property_ContigLocOutofInto :: Loc.ContigLoc -> Loc.ContigLoc -> Property
property_ContigLocOutofInto subcloc supercloc
    = let !mOutcloc = Loc.clocOutof subcloc supercloc
          !mIncloc = mOutcloc >>= flip Loc.clocInto supercloc
      in (isJust mOutcloc) ==> mIncloc == (Just subcloc)

test_Contig_IntoBounds :: Loc.ContigLoc -> Pos.Pos -> Bool
test_Contig_IntoBounds contig pos
    = let !mInpos = Loc.posInto pos contig
          !offset = Pos.offset pos
          !(cstart, cend) = Loc.bounds contig
      in (isJust mInpos) == (offset >= cstart && offset <= cend)

test_Contig_OutofBounds :: Loc.ContigLoc -> Pos.Pos -> Bool
test_Contig_OutofBounds contig pos
    = let !offset = Pos.offset pos
      in (isJust $ Loc.posOutof pos contig) == (offset >= 0 && offset < Loc.length contig)

property_Contig_allPos_outof :: Loc.ContigLoc -> Property
property_Contig_allPos_outof contig
  = forAll (choose (0, fromIntegral $ Loc.length contig - 1)) $ \ioff ->
  let p = drop ioff $ Loc.allPos contig
  in and [ not $ null p
         , Loc.posOutof (Pos.Pos (Pos.Offset $ fromIntegral ioff) Plus) contig == Just (head p)
         ]

property_Contig_seqData :: Loc.ContigLoc -> Property
property_Contig_seqData contig
    = forAll (genNonNegOffset >>= genNtByteString . fromIntegral) $ \sequ ->
      let seqData = Loc.seqData sequ contig
          padded = Loc.seqDataPad sequ contig
      in case seqData of
           (Just subsequ) -> and [ padded == subsequ, 'N' `BS.notElem` padded ]
           Nothing -> 'N' `BS.elem` padded

property_Contig_seqDataPadded :: Loc.ContigLoc -> Property
property_Contig_seqDataPadded contig
    = forAll (genNonNegOffset >>= genNtByteString . fromIntegral) $ \sequ ->
      (BS.pack $ map (fromMaybe 'N' . Pos.atPos sequ) contigPoses) == Loc.seqDataPad sequ contig
    where contigPoses = mapMaybe (flip Loc.posOutof contig . flip Pos.Pos Plus) [0..(Loc.length contig - 1)]

property_Contig_seqData2 :: Loc.ContigLoc -> Property
property_Contig_seqData2 contig
    = forAll (genNonNegOffset >>= genNtByteString . fromIntegral) $ \sequ ->
    let toLBS = LBS.fromChunks . (: [])
        fromLBS = BS.concat . LBS.toChunks
    in and [ Loc.seqData    sequ contig == liftM fromLBS (Loc.seqData    (toLBS     sequ) contig)
           , Loc.seqDataPad sequ contig ==       fromLBS (Loc.seqDataPad (toLBS     sequ) contig) 
           , Loc.seqData    sequ contig == liftM BS.pack (Loc.seqData    (BS.unpack sequ) contig)
           , Loc.seqDataPad sequ contig ==       BS.pack (Loc.seqDataPad (BS.unpack sequ) contig)
           ]

property_Contig_extendMinus :: Loc.ContigLoc -> Property
property_Contig_extendMinus contig
    = forAll (liftM2 (,) genNonNegOffset genNonNegOffset) $ \(ext5, ext3) ->
      (revCompl $ Loc.extend (ext5, ext3) contig) == (Loc.extend (ext3, ext5) $ revCompl contig)

property_Contig_fromStartEnd :: Loc.ContigLoc -> Property
property_Contig_fromStartEnd contig
    = (Loc.length contig > 1) ==>
      (Loc.fromStartEnd (Pos.offset $ Loc.startPos contig) (Pos.offset $ Loc.endPos contig)) == contig

property_Contig_fromBoundsStrand :: Loc.ContigLoc -> Property
property_Contig_fromBoundsStrand contig
    = (Loc.length contig > 1) ==>
      (Loc.fromBoundsStrand (fst . Loc.bounds $ contig) (snd . Loc.bounds $ contig) (Loc.strand contig)) == contig

property_Contig_overlaps :: Loc.ContigLoc -> Loc.ContigLoc -> Bool
property_Contig_overlaps cloc1 cloc2
  = (cloc1 `Loc.contigOverlaps` cloc2) ==
    and [ Loc.strand cloc1 == Loc.strand cloc2
        , or [ isJust . Loc.posInto (Loc.startPos cloc1) $ cloc2
             , isJust . Loc.posInto (Loc.startPos cloc2) $ cloc1
             , isJust . Loc.posInto (Loc.endPos   cloc1) $ cloc2
             , isJust . Loc.posInto (Loc.endPos   cloc2) $ cloc1
             ]
        ]

test_Contig_repr :: Loc.ContigLoc -> Bool
test_Contig_repr = test_repr

-- Bio.BioSeq.Location

genInvertibleLoc :: Gen SpLoc.SpliceLoc
genInvertibleLoc = sized $ \sz -> do ncontigs <- choose (1, sz + 1)
                                     fwdloc <- liftM (fromJust . SpLoc.fromContigs) 
                                               $ genContigs ncontigs
                                     rc <- arbitrary
                                     if rc then return $ revCompl fwdloc else return fwdloc
    where genContigs = liftM (reverse . foldl' intervalsToContigs []) . genIntervals
          genIntervals nints = replicateM nints $ liftM2 (,) genPositiveOffset genPositiveOffset
          intervalsToContigs [] (init5, len) = [Loc.fromPosLen (Pos.Pos init5 Plus) len]
          intervalsToContigs prevs@(prev:_) (nextoffset, nextlen)
              = let !prevend = Loc.offset5 prev + Loc.length prev
                in (Loc.fromPosLen (Pos.Pos (prevend + nextoffset) Plus) nextlen):prevs

instance Arbitrary SpLoc.SpliceLoc where
  arbitrary = genInvertibleLoc

test_Loc_Minus :: SpLoc.SpliceLoc -> Bool
test_Loc_Minus = test_revCompl

property_LocIntoOutof :: SpLoc.SpliceLoc -> Pos.Pos -> Property
property_LocIntoOutof loc pos
    = let !mInpos = Loc.posInto pos loc
          !mOutpos = mInpos >>= flip Loc.posOutof loc
      in (isJust mInpos) ==> mOutpos == (Just pos)

property_LocOutofInto :: Pos.Pos -> Property
property_LocOutofInto pos
    = forAll genInvertibleLoc $ \loc ->
      let !mOutpos = Loc.posOutof pos loc
          !mInpos = mOutpos >>= flip Loc.posInto loc
      in (isJust mOutpos) ==> mInpos == (Just pos)

test_Loc_OutofBounds :: SpLoc.SpliceLoc -> Pos.Pos -> Bool
test_Loc_OutofBounds loc pos
    = let !offset = Pos.offset pos
      in (isJust $ Loc.posOutof pos loc) == (offset >= 0 && offset < Loc.length loc)

property_LocCLocOutofInto :: Loc.ContigLoc -> Property
property_LocCLocOutofInto cloc
    = forAll genInvertibleLoc $ \loc ->
      let !mOutloc = Loc.clocOutof cloc loc
          !mInloc = mOutloc >>= mapM (flip Loc.clocInto loc) . Loc.toContigs >>= return . fromJust . SpLoc.fromContigs
      in (isJust mOutloc) ==> and [ liftM Loc.length mInloc == Just (Loc.length cloc)
                                  , liftM Loc.bounds mInloc == Just (Loc.bounds cloc)
                                  ]

property_LocOutofAssoc :: SpLoc.SpliceLoc -> Loc.ContigLoc -> Pos.Pos -> Property
property_LocOutofAssoc loc cloc pos
    = let !mOutloc = Loc.clocOutof cloc loc
          !mOutpos = mOutloc >>= \outloc -> Loc.posOutof pos outloc
      in (isJust mOutpos) ==> mOutpos == (Loc.posOutof pos cloc >>= \outpos -> Loc.posOutof outpos loc)

property_SpLocOutof :: SpLoc.SpliceLoc -> SpLoc.SpliceLoc -> Bool
property_SpLocOutof subloc outerloc =
  let !mOutofContigs = liftM Loc.toContigs $ SpLoc.locOutof subloc outerloc
      !mContigsOutof = liftM (concat . map Loc.toContigs) $
                       mapM (flip Loc.clocOutof outerloc) $
                       Loc.toContigs subloc
  in mOutofContigs == mContigsOutof

property_SpLocOutofGood :: SpLoc.SpliceLoc -> SpLoc.SpliceLoc -> Property
property_SpLocOutofGood subloc outerloc =
  let !mOutofContigs = liftM Loc.toContigs $ SpLoc.locOutof subloc outerloc
      !mContigsOutof = liftM (concat . map Loc.toContigs) $
                       mapM (flip Loc.clocOutof outerloc) $
                       Loc.toContigs subloc
  in (isJust mOutofContigs) ==> mOutofContigs == mContigsOutof

property_Loc_allPos_outof :: SpLoc.SpliceLoc -> Property
property_Loc_allPos_outof sploc
  = forAll (choose (0, fromIntegral $ Loc.length sploc - 1)) $ \ioff ->
  let p = drop ioff $ Loc.allPos sploc
  in and [ not $ null p
         , Loc.posOutof (Pos.Pos (Pos.Offset $ fromIntegral ioff) Plus) sploc == Just (head p)
         ]

property_Loc_seqData :: SpLoc.SpliceLoc -> Property
property_Loc_seqData loc
    = forAll (genNonNegOffset >>= genNtByteString . fromIntegral) $ \sequ ->
      let seqData = Loc.seqData sequ loc
          padded = Loc.seqDataPad sequ loc
      in case seqData of
           (Just subsequ) -> and [ padded == subsequ, 'N' `BS.notElem` padded ]
           Nothing -> 'N' `BS.elem` padded

property_Loc_seqDataPadded :: SpLoc.SpliceLoc -> Property
property_Loc_seqDataPadded loc
    = forAll (genNonNegOffset >>= genNtByteString . fromIntegral) $ \sequ ->
      (BS.pack $ map (fromMaybe 'N' . Pos.atPos sequ) locPoses) == Loc.seqDataPad sequ loc
    where locPoses = mapMaybe (flip Loc.posOutof loc . flip Pos.Pos Plus) [0..(Loc.length loc - 1)]

property_SpLoc_seqData2 :: SpLoc.SpliceLoc -> Property
property_SpLoc_seqData2 sploc
    = forAll (genNonNegOffset >>= genNtByteString . fromIntegral) $ \sequ ->
    let toLBS = LBS.fromChunks . (: [])
        fromLBS = BS.concat . LBS.toChunks
    in and [ Loc.seqData    sequ sploc == liftM fromLBS (Loc.seqData    (toLBS     sequ) sploc)
           , Loc.seqDataPad sequ sploc ==       fromLBS (Loc.seqDataPad (toLBS     sequ) sploc) 
           , Loc.seqData    sequ sploc == liftM BS.pack (Loc.seqData    (BS.unpack sequ) sploc)
           , Loc.seqDataPad sequ sploc ==       BS.pack (Loc.seqDataPad (BS.unpack sequ) sploc)
           ]

property_Loc_Within :: Pos.Pos -> Property
property_Loc_Within pos
    = forAll genInvertibleLoc $ \loc ->
      and [ (pos `Loc.posWithin` loc) == (maybe False ((/= Minus) . Pos.strand) $ Loc.posInto pos loc)
          , ((Pos.offset pos) `Loc.offsetWithin` loc) == (isJust . Loc.posInto pos $ loc)
          ]

test_SpLoc_repr :: SpLoc.SpliceLoc -> Bool
test_SpLoc_repr = test_repr

property_SpLoc_terminiMinus :: SpLoc.SpliceLoc -> Bool
property_SpLoc_terminiMinus loc
  = and [ revCompl (Loc.startPos loc) == Loc.endPos (revCompl loc)
        , revCompl (Loc.endPos loc) == Loc.startPos (revCompl loc)
        ]

property_SpLoc_terminiExtend :: SpLoc.SpliceLoc -> Property
property_SpLoc_terminiExtend loc  
  = forAll genNonNegOffset $ \ext5 ->
  forAll genNonNegOffset $ \ext3 ->
  let extloc = Loc.extend (ext5, ext3) loc
      strandSlide pos doff = case Pos.strand pos of
        Plus -> Pos.slide pos doff
        Minus -> Pos.slide pos (negate doff)
  in and [ Loc.startPos extloc == strandSlide (Loc.startPos loc) (negate ext5)
         , Loc.endPos extloc == strandSlide (Loc.endPos loc) ext3
         ]

property_ShVector_singleton :: Property
property_ShVector_singleton =
  forAll (liftM fromIntegral genOffset) $ \i ->
  forAll (liftM fromIntegral genOffset) $ \j ->
  forAll arbitrary $ \(ch :: Char) ->
  let sv = ShV.singleton i [ch]
  in and [ sv ShV.!? i == [ch],
           sv ShV.!? (i - 1) == [],
           sv ShV.!? (i + 1) == [],
           sv ShV.!? j == if (i == j) then [ch] else [] ]

property_ShVector_update1 :: Property
property_ShVector_update1 = 
  forAll (liftM fromIntegral genOffset) $ \i ->
  forAll arbitrary $ \(chi :: Char) ->
  let sv0 = ShV.empty
      sv1 = sv0 ShV.// [(i, [chi])]
  in and [ sv1 ShV.!? i == [chi],
           sv1 ShV.!? (i - 1) == [],
           sv1 ShV.!? (i + 1) == [] ]

property_ShVector_update2 :: Property
property_ShVector_update2 = 
  forAll (liftM fromIntegral genOffset) $ \i ->
  forAll (liftM fromIntegral genOffset) $ \j ->
  forAll arbitrary $ \(chi :: Char) ->
  forAll arbitrary $ \(chj :: Char) ->
  let sv0 = ShV.empty
      sv1 = sv0 ShV.// [(i, [chi]), (j, [chj])]
  in (i /= j) ==>
     and [ sv1 ShV.!? i == [chi],
           sv1 ShV.!? j == [chj],
           sv1 ShV.!? ((min i j) - 1) == [],
           sv1 ShV.!? ((max i j) + 1) == [],
           and [ sv1 ShV.!? k == [] | k <- [(min i j + 1)..(max i j - 1)] ] ]

property_LocMap_hitIn :: Loc.ContigLoc -> Pos.Pos -> Property
property_LocMap_hitIn contig pos =
  forAll genPositiveOffset $ \binsz -> 
  let isin = isJust $ Loc.posInto pos contig
      ploc = Loc.fromPosLen pos 1
      lm = LM.insertLoc contig contig (LM.emptyLM binsz)
  in collect isin $ isin ==> not (null (LM.queryLoc ploc lm))

property_LocMap_hitAll :: Loc.ContigLoc -> Pos.Pos -> Property
property_LocMap_hitAll contig pos =
  forAll genPositiveOffset $ \binsz -> 
  let isin = isJust $ Loc.posInto pos contig
      ploc = Loc.fromPosLen pos 1
      lm = LM.insertLoc contig contig (LM.emptyLM binsz)
      ishit = not $ null (LM.queryLoc ploc lm)
  in collect ishit $ ishit || not isin

property_LocMap_hitMulti :: Loc.ContigLoc -> Loc.ContigLoc -> Loc.ContigLoc -> Pos.Pos -> Property
property_LocMap_hitMulti ca cb cc pos =
  forAll genPositiveOffset $ \binsz -> 
  let isina = isJust $ Loc.posInto pos ca
      isinb = isJust $ Loc.posInto pos cb
      isinc = isJust $ Loc.posInto pos cc
      ploc = Loc.fromPosLen pos 1
      lm = LM.insertLoc cc cc $
           LM.insertLoc cb cb $
           LM.insertLoc ca ca (LM.emptyLM binsz)
      hita = ca `elem` LM.queryLoc ploc lm
      hitb = cb `elem` LM.queryLoc ploc lm
      hitc = cc `elem` LM.queryLoc ploc lm
  in and [ hita || not isina,
           hitb || not isinb,
           hitc || not isinc ]

property_Locatable_hitWithin :: Loc.ContigLoc -> Loc.ContigLoc -> Pos.Pos -> Property
property_Locatable_hitWithin ca cb pos =
  forAll genPositiveOffset $ \binsz ->
  forAll genName $ \chr ->
  let isina = isJust $ Loc.posInto pos ca
      isinb = isJust $ Loc.posInto pos cb
      ploc = OnSeq chr $ Loc.fromPosLen pos 1
      la = LM.WithLocation ca (OnSeq chr (SpLoc.singleton ca))
      lb = LM.WithLocation cb (OnSeq chr (SpLoc.singleton cb))
      lm = LM.locatableSeqLocMap binsz [ la, lb ]
      hita = la `elem` LM.queryLocatable Nothing ploc lm
      hitb = lb `elem` LM.queryLocatable Nothing ploc lm
      hita' = la `elem` LM.querySeqLoc ploc lm
      hitb' = lb `elem` LM.querySeqLoc ploc lm
  in and [ hita == isina, hita' || not isina
         , hitb == isinb, hitb' || not isinb
         ]

instance (Arbitrary a) => Arbitrary (OnSeq a) where
  arbitrary = OnSeq <$> arbitrary <*> arbitrary

-- Utilities

data Test = forall t . Testable t => T String t

runTest :: Test -> IO ()
runTest (T name test) = do
  putStr $ name ++ replicate (40 - length name) '.' ++ "  "
  quickCheckWith args test
    where args = stdArgs -- { maxDiscard = 100000 }

-- | Constrained position generators

genOffset :: Gen Pos.Offset
genOffset = do isneg <- arbitrary
               nnoff <- genNonNegOffset
               return $ (if isneg then negate else id) nnoff

genNonNegOffset :: Gen Pos.Offset
genNonNegOffset = liftM (subtract 1) genPositiveOffset

genPositiveOffset :: Gen Pos.Offset
genPositiveOffset = do scale <- chooseInteger (1, 10)
                       liftM fromIntegral $ chooseInteger (1, 2^scale)
    where chooseInteger :: (Integer, Integer) -> Gen Integer
          chooseInteger = choose

