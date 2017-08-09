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

## Requirements

To run a pipeline, please provide a configuration file with the paths to the various tool directories and executables.  Please see the [example configuration file](https://github.com/fulcrumgenomics/fgstr/blob/devel/pipelines/src/main/resources/application.conf).
More generally, the following tools are required:

  - [fgbio](https://github.com/fulcrumgenomics/fgbio)
  - [fgstr](https://github.com/fulcrumgenomics/fgstr)
  - [picard](https://github.com/broadinstitute/picard/)
  - [bwa](https://github.com/lh3/bwa)
  - [HipSTR](https://github.com/tfwillems/HipSTR)

## STR Calling a Single Sample with Duplex Sequencing

### Quick Start

The following command should be used for single sample analysis:

`java -jar pipelines/target/scala-2.12/fgstr-pipelines-0.1.0-SNAPSHOT.jar GenotypeWithDuplexSequencing ...` 

## Background

STR calling is performed by treating the reads from each duplex source molecule as a single haploid sample.
The set of most frequent alleles are called as the genotype, such that the total frequency is greater than a fixed 
threshold (ex. `0.9`)

To analyze a single sample use the `GenotypeWithDuplexSequencing` pipeline.  It automates the following three pipelines:
  1. `PrepareUnmappedBamPipeline`: create an unmapped BAM with extracted UMI information using various [picard](https://github.com/broadinstitute/picard/)
  and [fgbio](https://github.com/fulcrumgenomics/fgbio) tools.
  2. `MapAndGroupRawReads`: map the raw reads (with [bwa](https://github.com/lh3/bwa) and assign them each to a single source duplex molecule (with
  [fgbio](https://github.com/fulcrumgenomics/fgbio)'s `GroupReadsByUmi`). Additionally, collects various QC metrics.
  3. `GenotypeFromGroupedBam`: generate STR genotypes from the grouped BAM.  This entails a) generating a BAM with a 
  read-group per source molecule, b) running [HipSTR](https://github.com/tfwillems/HipSTR) to jointly call across all
  duplex source molecules (treating each source as a single haploid sample), and c) combining the per-molecule "genotypes" 
  to produce a final genotype.

The details of each step can be found in the pipeline sources, and key tools to review are:
  - [fgbio](https://github.com/fulcrumgenomics/fgbio)'s `GroupReadsByUmi`
  - [fgstr](https://github.com/fulcrumgenomics/fgstr)'s `CollectDuplexSeqMetrics`
  - [fgstr](https://github.com/fulcrumgenomics/fgstr)'s `ReadGroupPerDuplexMolecularId`
  - [HipSTR](https://github.com/tfwillems/HipSTR) main executable
  - [HipSTR](https://github.com/tfwillems/HipSTR)'s `filter_haploid_vcf.py` script
  - [fgstr](https://github.com/fulcrumgenomics/fgstr)'s `StrGenotypeDuplexMolecules`

### Input Requirements

The interval list provided to both `GenotypeWithDuplexSequencing ` and `GenotypeFromGroupedBam` requires a custom format.
The name field should contain a comma list of values as follows:

  1. the reference genome repeat unit length (ex. `3` for the tri-nucleotide repeat `TCATCATCATCA`).
  2. the reference genome number of repeat units (ex. `4` for the tri-nucleotide repeat `TCATCATCATCA`).
  3. the name of the STR (ex. `D1S1656`)
 
Optionally, additional values can be appended to the name for one or more expected truth alleles.  
For example, a known haploid call should have one extra value, a known diploid call should have two extra value, and so on.
  
An example name field with the optionals is `4,17,D1S1656,9,10` and without is `4,17,D1S1656`.

For example, below is an interval list with (1) the header (sequence dictionary) omitted, and dummy values for the expected/known/truth number of repeat units.  
Notice that `D12S391` is tri-allelic.

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
12  12449955    12450029    +   4,19,D12S391,9,10,11
13  82722161    82722203    +   4,11,D13S317,9,10
15  97374246    97374269    +   6,5,PentaE,9,10
16  86386309    86386351    +   4,11,D16S539,9,10
18  60948901    60948971    +   4,18,D18S51,9,10
19  30417143    30417205    +   4,14,D19S433,9,10
21  20554292    20554417    +   4,29,D21S11,29,31.2
21  43636206    43636269    +   5,13,PentaD,9,10
22  37536328    37536377    +   3,17,D22S1045,9,10
```

### Multi-sample processing

Multiple samples can be processed using the`MultiSampleGenotypeWithDuplexSequencing` pipeline. 
The following directory structure should be present as input:
```
root
|--inputs
|  |--sample.1
|     |--..._L001_R1_001.fastq.gz
|     |--..._L001_R2_001.fastq.gz
|  |...
|  |--sample.n
|     |--..._L001_R1_001.fastq.gz
|     |--..._L001_R2_001.fastq.gz
|
|--interval_lists
   |--sample.1.interval_list
   |...
   |--sample.n.interval_list
```
