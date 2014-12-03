Introduction
============

These auxilliary `Bio.SeqLoc` modules provide functions to read and
write gene annotations in GTF and BED format.  I separated them to
minimize the number of dependencies for the main `seqloc` package.

GTF
===

A single GTF annotation can span multiple lines and the "[o]rder of
rows is not important", so the entire file must be loaded before any
transcripts can be assembled.

The following program reads a full human GTF annotation, takes the
first 10 transcripts, writes them to a small test file, re-reads them,
and verifies that the results of a cycle of `transcriptToGtf` followed
by `readGtfTranscripts` does not change anything.

	module Main
	       where
	
	import Control.Monad
	import qualified Data.ByteString.Char8 as BS
	import Data.List
	import System.IO
	
	import Bio.SeqLoc.GTF
	import Bio.SeqLoc.LocRepr
	import Bio.SeqLoc.OnSeq
	import Bio.SeqLoc.Transcript
	
	main :: IO ()
	main = do trx <- readGtfTranscripts "/data/genomes/Homo_sapiens/hg19_knownGene.gtf"
	          let trx10 = take 10 trx'
	          BS.writeFile "test/gtf-out10.gtf" . BS.concat . map (transcriptToGtf "TestGtf") $ trx10
	          trx10' <- readGtfTranscripts "test/gtf-out10.gtf"
	          print $ (sort . map location $ trx10) == (sort . map location $ trx10')

BED
===

A single BED annotation occupies a single line, so it is possible to
process BED format annotations iteratively. This interface uses
`Data.Iteratee` iteration, as shown below.

The `Iter.mapM_` function generates an `Iteratee` that maps a monadic
action over each element of the input stream. Here, the input stream
will be a list of transcripts, which will be written to the output
file using `BS.hPutStrLn hout . transcriptToBedStd`. The
`bedTranscriptEnum` encloses the BED format parser, allowing it to
convert an iteratee for a stream of `[Transcript]` into an iteratee
for a stream of `BS.ByteString` containing the BED format data. The
`Iter.fileDriver` driver then applies the transformed `bedIter` to the
data from a file.

	module Main
	       where
	
	import Control.Monad
	import qualified Data.ByteString.Char8 as BS
	import Data.List
	import System.IO
	
	import qualified Data.Iteratee as Iter
	
	import Bio.SeqLoc.Bed
	import Bio.SeqLoc.LocRepr
	import Bio.SeqLoc.OnSeq
	import Bio.SeqLoc.Transcript
	
	main :: IO ()
	main = do withFile "test/bed-copy.bed" WriteMode $ \hout ->
	            let bedIter = bedTranscriptEnum $ Iter.mapM_ (BS.hPutStrLn hout . transcriptToBedStd)
	            in Iter.fileDriver bedIter "/data/genomes/Homo_sapiens/hg19_knownGene.bed"

A simpler interface in which the entire contents of an annotation file
are read or written together is also provided, just as for GTF.
