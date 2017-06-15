[![Build Status](https://travis-ci.org/fulcrumgenomics/fgstr.svg?branch=devel)](https://travis-ci.org/fulcrumgenomics/fgstr)
[![codecov](https://codecov.io/gh/fulcrumgenomics/fgstr/branch/devel/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fgstr)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/ff514740e09f4ed0a62fbff805a2de59)](https://www.codacy.com/app/nilshomer/fgstr?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=fulcrumgenomics/fgstr&amp;utm_campaign=Badge_Grade)[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fgstr/blob/devel/LICENSE)
[![Language](http://img.shields.io/badge/language-scala-brightgreen.svg)](http://www.scala-lang.org/)

# fgstr

Simple Tandem Repeat Genotyping for Duplex Sequencing

## *** WARNING ***

This tools and pipelines in this repository are in active development, and are not intended for production use.
They are subject to change at any time as we develop new methods and procedures.

## Building 
### Cloning the Repository

`git clone https://github.com/fulcrumgenomics/fgstr.git`

### Running the Build
fgstr is built using [sbt](http://www.scala-sbt.org/).

Use ```sbt assembly``` to build an executable jar in ```target/scala-2.12/```.  
Tests may be run with ```sbt test```.
Java SE 8 is required.

## STR Calling a Single Sample with Duplex Sequencing

STR calling is performed by treating the reads from each duplex source molecule as a single haploid sample.  
The two most frequently called set of alleles are called as the genotype for the sample.

To analyze a single sample:

1. Run `PrepareUnmappedBamPipeline` to create an unmapped BAM with extracted UMI information.
2. Run `MapAndGroupRawReads` to map the raw reads and assign them each to a single source duplex molecule.
3. Run `GenotypeFromGroupedBam` to generate STR genotypes from the grouped BAM. 
This entails a) generating a BAM with a read-group per source molecule, b) running HipSTR to jointly call across all duplex source molecules (treating each source as a single haploid sample), and c) combining the per-molecule "genotypes" to produce a final genotype.

### Input Requirements

The interval list provided to `GenotypeFromGroupedBam` requires a custom format.
The name field should contain a comma list of values as follows:

  1. the reference genome repeat unit length (ex. `3` for the tri-nucleotide repeat `TCATCATCATCA`).
  2. the reference genome number of repeat units (ex. `4` for the tri-nucleotide repeat `TCATCATCATCA`).
  3. the name of the STR (ex. `D1S1656`)
 
Optionally, two extra values can be appended:

  4. optionally, the expected (known or truth) number of repeat units for allele #1
  5. optionally, the expected (known or truth) number of repeat units for allele #2
  
An example name field with the optionals is `4,17,D1S1656,9,10` and without is `4,17,D1S1656`.

For example, below is an interval list with (1) the header (sequence dictionary) omitted, and (2) dummy values for the expected/known/truth number of repeat units:

```
1   230905363   230905429   +   4,17,D1S1656,9,10
2   1493426 1493456 +   4,8,TPOX,9,10
2   68239080    68239126    +   4,12,D2S441,9,10
2   218879583   218879673   +   4,23,D2S1338,9,10
3   45582232    45582294    +   4,16,D3S1358,9,10
4   155508889   155508975   +   4,22,FGA,9,10
5   123111251   123111293   +   4,11,D5S818,9,10
5   149455888   149455938   +   4,13,CSF1PO,9,10
6   88986849    88986964    +   4,25,SE33,9,10
6   92449944    92449990    +   4,12,D6S1043,9,10
7   83789543    83789593    +   4,13,D7S820,9,10
8   125907108   125907158   +   4,13,D8S1179,9,10
10  131092509   131092559   +   4,13,D10S1248,9,10
11  2192319 2192345 +   4,7,TH01,6,9.3
12  6093144 6093210 +   4,17,vWA,9,10
12  12449955    12450029    +   4,19,D12S391,9,10
13  82722161    82722203    +   4,11,D13S317,9,10
15  97374246    97374269    +   6,5,PentaE,9,10
16  86386309    86386351    +   4,11,D16S539,9,10
18  60948901    60948971    +   4,18,D18S51,9,10
19  30417143    30417205    +   4,14,D19S433,9,10
21  20554292    20554417    +   4,29,D21S11,29,31.2
21  43636206    43636269    +   5,13,PentaD,9,10
22  37536328    37536377    +   3,17,D22S1045,9,10
```
