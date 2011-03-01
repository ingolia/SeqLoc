module Bio.SeqLoc.Transcript
       (
         -- * Type for splice junctions
         Junction (..)
       , fromDonorAcceptor, donor, acceptor
       , junctions
         -- * Representation of transcript
       , Transcript(..)         
       , cdsLocation
       , sortContigs
       )
       where 

import Control.Applicative
import Control.Monad
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Ord

import qualified Data.Attoparsec.Zepto as ZP

import Bio.SeqLoc.LocRepr
import qualified Bio.SeqLoc.Location as Loc
import Bio.SeqLoc.OnSeq
import qualified Bio.SeqLoc.Position as Pos
import qualified Bio.SeqLoc.SpliceLocation as SpLoc
import Bio.SeqLoc.Strand

-- | Splice junctions, which are isomorphic to the introns they span,
-- but which support other biologically relevant constructors and
-- accessors.
newtype Junction = Junction { intron :: Loc.ContigLoc } deriving (Show)

slash :: BS.ByteString
slash = BS.pack "/"

instance LocRepr Junction where
  repr j = BS.concat [ repr . donor $ j, slash, repr . acceptor $ j ]
  unrepr = fromDonorAcceptor <$> unrepr <*> (ZP.string slash *> unrepr)  

-- | Create a splice junction from a donor position (the last position
-- in the 5' exon) and the acceptor position (the first position in
-- the 3' exon).
fromDonorAcceptor :: Pos.Pos -> Pos.Pos -> Junction
fromDonorAcceptor d a = let len = 1 + abs (Pos.offset a - Pos.offset d)
                        in case Pos.strand d of
                          Fwd -> Junction $! Loc.fromPosLen (Pos.slide d 1) len
                          RevCompl -> Junction $! Loc.fromPosLen (Pos.slide d (-1)) len

-- | Donor position, i.e., the last position in the 5' exon around a
-- junction.
donor :: Junction -> Pos.Pos
donor = Loc.startPos . Loc.extend (1, 0) . intron

-- | Acceptor position, i.e., the first position in the 3' exon around
-- a junction.
acceptor :: Junction -> Pos.Pos
acceptor = Loc.endPos . Loc.extend (0, 1) . intron

-- | List of splice junctions from a spliced location, in order.
junctions :: SpLoc.SpliceLoc -> [Junction]
junctions sploc = zipWith junction contigs (drop 1 contigs)
  where contigs = Loc.toContigs sploc
        junction c5 c3 = let p5 = Loc.endPos . Loc.extend (0, 1) $ c5
                             p3 = Loc.startPos . Loc.extend (1, 0) $ c3
                             len = 1 + abs (Pos.offset p3 - Pos.offset p5)
                         in Junction $ Loc.fromPosLen p5 len

-- | Representation of a genomic transcript, with a gene and a
-- transcript identifier, along with the genomic location of the
-- processed transcript and an optional coding sequence on that
-- transcript.
data Transcript = Transcript { geneId -- ^ Gene or locus name for a collection of transcripts
                             , trxId :: !SeqName -- ^ Specific transcript identifier
                             , location :: !SpliceSeqLoc -- ^ Sequence location of processed transcript
                             , cds :: !(Maybe Loc.ContigLoc) -- ^ Location of CDS on the transcript
                             } deriving (Show)
                                        
-- | Genomic location of CDS within the transcript
cdsLocation :: Transcript -> Maybe SpliceSeqLoc
cdsLocation trx = cds trx >>= liftM (OnSeq name) . flip Loc.clocOutof loc
  where (OnSeq name loc) = location trx

-- | 'Just' the input contigs sorted in stranded order, when all lie
-- on the same strand, or 'Nothing' if they are not all on the same
-- strand.
sortContigs :: [Loc.ContigLoc] -> Maybe [Loc.ContigLoc]
sortContigs [] = Nothing
sortContigs cs@(c0:_)= liftM sortStrand contigStrand
  where contigStrand | all ((== Loc.strand c0) . Loc.strand) cs = Just . Loc.strand $ c0
                     | otherwise = Nothing
        sortStrand Fwd = sortBy (comparing Loc.offset5) cs
        sortStrand RevCompl = sortBy (comparing (negate . Loc.offset5)) cs