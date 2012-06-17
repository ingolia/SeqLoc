Introduction
============

The `Bio.SeqLoc` modules in `seqloc` are designed to represent
positions and locations (ranges of positions) on sequences,
particularly nucleotide sequences. My original motivation for writing
these packages was handing the locations of genes in eukaryotic
genomes.

Strands
=======

Handling forward and reverse-complement locations and sequences is a
very common task in bioinformatics. The `Bio.SeqLoc.Strand` package
handles strandedness in `seqloc`. It consists of a simple `Strand`
enumerated data type and a `Stranded` typeclass of any object that has
a strandedness and can therefore be meaningfully reverse complemented.
String-like objects such as `String` itself and the various
`ByteString` types instantiate `Stranded` by reversing the string and
complementing nucleotide characters.

	Prelude> :set prompt "ghci> "
	ghci> :set -XOverloadedStrings
	ghci> import Bio.SeqLoc.Strand
	ghci> revCompl "GAttACA"
	"TGTaaTC"
	ghci> stranded Plus "GATTaca"
	"GATTaca"
	ghci> stranded Minus "GATTaca"
	"tgtAATC"
	ghci> stranded (revCompl Plus) "GATTaca"
	"tgtAATC"

Positions and Offsets
=====================

	ghci> import Bio.SeqLoc.Position

In the `seqloc` package, an `Offset` is a 0-based index into a
sequence and a `Position` is an `Offset` plus a `Strand` indicating
the strand on which the position occurs in the sequence.

Displaying Positions and Locations
==================================

	ghci> import Bio.SeqLoc.LocRepr
	ghci> :m +Data.ByteString.Char8

The `LocRepr` typeclass provides an interface for representing
position and location data types in a format that is easy to read as
well as to parse.  There are two basic functions, `repr` which
produces a string representation and `unrepr` which is an extremely
lightweight parser for that string representation from
`Data.Attoparsec.Zepto`.  There are also helper functions that wrap
the parser and handle errors in different ways.

	ghci> repr (Pos 99 Plus)
	"99(+)"
	ghci> (unreprEither  "99(-)") :: Either String Pos
	Right (Pos {offset = Offset {unOffset = 99}, strand = Minus})
	ghci> (unreprErr  "99(+)") :: Pos
	Pos {offset = Offset {unOffset = 99}, strand = Plus}


While locations and positions have `Show` instances, their `LocRepr`
instances have advantages for human and computer legibility in many
contexts.

Contiguous Locations
====================

	ghci> import Bio.SeqLoc.Location as Loc

The `ContigLoc` type represents a contiguous sequence location, such
as the forward strand from nucleotides 100 to 150, or the reverse
complement strand from nucleotides 1000 to 800. These locations can be
created by specifying their bounds and strand.

	ghci> repr $ Loc.fromBoundsStrand 100 150 Plus
	"100to150(+)"
	ghci> repr $ Loc.fromBoundsStrand 800 1000 Minus
	"800to1000(-)"

They can also be specified with a starting position, which will be the
beginning of the location *in its strand*, and a length.

	ghci> repr $ Loc.fromPosLen (unreprErr "100(+)") 51
	"100to150(+)"
	ghci> repr $ Loc.fromPosLen (unreprErr "1000(-)") 201
	"800to1000(-)"

Finally, they can be specified from their starting and ending
position, in which case the strand is deduced from the order of the
two positions.

	ghci> repr $ Loc.fromStartEnd 100 150
	"100to150(+)"
	ghci> repr $ Loc.fromStartEnd 1000 800
	"800to1000(-)"

The `ContigLoc` type is an instance of the `Location` typeclass, which
provides numerous useful functions.  It's important to remember that
the starting position of a (-) strand location has a higher offset
than the ending position.

	ghci> let l = Loc.fromStartEnd 1000 800
	ghci> let l' = revCompl l
	ghci> Loc.bounds l
	(Offset {unOffset = 800},Offset {unOffset = 1000})
	ghci> repr $ Loc.startPos l
	"1000(-)"
	ghci> repr $ Loc.endPos l
	"800(-)"
	ghci> Loc.bounds l'
	(Offset {unOffset = 800},Offset {unOffset = 1000})
	ghci> repr $ Loc.startPos l'
	"800(+)"
	ghci> repr $ Loc.endPos l'
	"1000(+)"

The `Location` typeclass also allows us to convert between a position
in absolute coordinates and a position relative to a location. The
`posInto` function takes an absolute position *into* a
location-relative position. It may fail, with `Nothing`, if the
position is outside the location.

	ghci> let p = Pos 850 Plus
	ghci> maybe "n/a" repr $ Loc.posInto p l
	"150(-)"
	ghci> maybe "n/a" repr $ Loc.posInto p l'
	"50(+)"
	ghci> let p2 = Pos 750 Plus
	ghci> maybe "n/a" repr $ Loc.posInto p2 l
	"n/a"

The `posOutof` function pulls a location-relative position back out of
the location.

	ghci> let q = Pos 120 Plus
	ghci> maybe "n/a" repr $ Loc.posOutof q l
	"880(-)"
	ghci> maybe "n/a" repr $ Loc.posOutof q l'
	"920(+)"
	ghci> let q2 = Pos 220 Plus
	ghci> maybe "n/a" repr $ Loc.posOutof q2 l'
	"n/a"

Entire locations can be mapped back in a similar way. Here we find a
sub-location from nucleotides 100 through 150 within the enclosing
location of 1000 to 800, map the sub-location back to its absolute
coordinates, and then find its relative coordinates within the
complementary location.

	ghci> let k = Loc.fromStartEnd 100 150
	ghci> maybe "n/a" repr $ Loc.clocOutof k l
	"850to900(-)"
	ghci> maybe "n/a" repr $ Loc.clocInto (unreprErr "850to900(-)") l'
	"50to100(-)"

Sequence Data
=============

	ghci> import Bio.SeqLoc.SeqLike as SeqLike

The `SeqLike` typeclass in the `Bio.SeqLoc.SeqLike` module has a
simple interface to allow the extraction of subsequences based on
locations. There are instances for `String` and for lazy and strict
`ByteString` types. Recall that offsets are all 0-based indices and
that `fromStartEnd` creates a location that includes both endpoints.

	ghci> Loc.seqData "GATTACA" (Loc.fromStartEnd 2 4)
	Just "TTA"
	ghci> Loc.seqData "GATTACA" (Loc.fromStartEnd 2 8)
	Nothing
	ghci> Loc.seqDataPad "GATTACA" (Loc.fromStartEnd 2 8)
	"TTACANN"
	ghci> Loc.seqDataPad "GATTACA" (Loc.fromStartEnd (-2) 4)
	"NNGATTA"
	ghci> Loc.seqDataPad "GATTACA" (Loc.fromStartEnd 6 0)
	"TGTAATC"

The instances for `String` and lazy `ByteString` avoid evaluating the
full sequence whenever possible, but the use of functions such as
`length` will force its evaluation.

Spliced Locations
=================

	ghci> import Bio.SeqLoc.SpliceLocation as SpLoc

The `SpliceLoc` type in the `Bio.SeqLoc.SpliceLocation` package
provides spliced locations, designed to model the structure of
eukaryotic genes as a series of individual `ContigLoc` locations lying
in order on the same strand.

	ghci> maybe "n/a" repr $ SpLoc.fromContigs [ Loc.fromStartEnd 100 150, Loc.fromStartEnd 200 250 ]
	"100to150(+);200to250(+)"
	ghci> maybe "n/a" repr $ SpLoc.fromContigs [ Loc.fromStartEnd 100 150, Loc.fromStartEnd 250 200 ]
	"n/a"
	ghci> maybe "n/a" repr $ SpLoc.fromContigs [ Loc.fromStartEnd 100 150, Loc.fromStartEnd 50 80 ]
	"n/a"

`SpliceLoc` implements the `Location` interface as well. When pulling
a contiguous location out of a spliced location, the result may also
be spliced. When pushing a contiguous location into a spliced
location, it must fit entirely within a single segment of the spliced
location.

	ghci> let (Just s) = SpLoc.fromContigs [ Loc.fromStartEnd 100 150, Loc.fromStartEnd 200 250 ]
	ghci> maybe "n/a" repr $ Loc.clocOutof (Loc.fromStartEnd 25 75) s
	"125to150(+);200to224(+)"
	ghci> maybe "n/a" repr $ Loc.clocInto (Loc.fromStartEnd 210 240) s
	"61to91(+)"
	ghci> maybe "n/a" repr $ Loc.clocInto (Loc.fromStartEnd 190 240) s
	"n/a"

The module also provides specialized functions for finding the
coordinates of spliced locations relative to an enclosing spliced
location. These will properly merge locations whose relative
coordinates are adjacent

	ghci> let (Just u) = Loc.clocOutof (Loc.fromStartEnd 20 81) s
	ghci> repr u
	"120to150(+);200to230(+)"
	ghci> let (Just t) = SpLoc.fromContigs [ Loc.fromStartEnd 120 150, Loc.fromStartEnd 200 230 ]
	ghci> maybe "n/a" repr $ SpLoc.locInto t s
	"20to81(+)"
	ghci> let (Just t') = SpLoc.fromContigs [ Loc.fromStartEnd 120 150, Loc.fromStartEnd 210 230 ]
	ghci> maybe "n/a" repr $ SpLoc.locInto t' s
	"20to50(+);61to81(+)"

Named Sequences
===============

	ghci> import Bio.SeqLoc.OnSeq
	ghci> import qualified Data.ByteString.Char8 as BS

Genome annotation data files typically express the location of a gene
as a [spliced] location on one of several chromosomes. The `OnSeq`
type in the `Bio.SeqLoc.OnSeq` allows position and location types to
be tagged with names. The module provides useful type synonyms for the
named location data types.

	ghci> let z = OnSeq (toSeqLabel "chr1") (Loc.fromStartEnd 10000 20000)
	ghci> repr z
	"chr1@10000to20000(+)"
	ghci> repr $ revCompl z
	"chr1@10000to20000(-)"

