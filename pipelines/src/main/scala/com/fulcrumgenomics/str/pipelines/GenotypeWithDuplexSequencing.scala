/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package com.fulcrumgenomics.str.pipelines

import com.fulcrumgenomics.commons.CommonsDef.DirPath
import com.fulcrumgenomics.sopt.{arg, clp}
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.Pipeline
import dagr.tasks.DagrDef.{PathPrefix, PathToFasta, PathToFastq, PathToIntervals}

@clp(group=classOf[Pipelines], description=
  """
    |Pipeline to genotype STRs using Duplex Sequencing.
    |
    |Automates the following pipelines:
    |  1. `PrepareUnmappedBamPipeline`
    |  2. `MapAndGroupRawReads`
    |  3. `GenotypeFromGroupedBam`
    |
    |## Inputs
    |
    |An interval list specifying the set of regions over which to call STRs should be given. The name field should
    |contain a comma list of values as follows:
    |  1. the repeat unit length (ex. `3` for the tri-nucleotide repeat `TCATCATCATCA`).
    |  2. the number of repeat units (ex. `4` for the tri-nucleotide repeat `TCATCATCATCA`).
    |  3. the name of the STR (ex. D1S1656)
    |  4. optionally, the expected (known or truth) number of repeat units for allele #1
    |  5. optionally, the expected (known or truth) number of repeat units for allele #2
    |An example name field with the optionals is `4,17,D1S1656,9,10` and without is `4,17,D1S1656`.
    |See [the wiki](https://github.com/fulcrumgenomics/fgstr#input-requirements) for more details.
    |
    |Please see the respective pipelines for more information on specific command line options and behavior.
    """)
class GenotypeWithDuplexSequencing
(
  @arg(flag='1', doc="Input fastq file (optionally gzipped) for read 1.") val fastq1: List[PathToFastq],
  @arg(flag='2', doc="Input fastq file (optionally gzipped) for read 2.") val fastq2: List[PathToFastq],
  @arg(flag='o', doc="Output file prefix (e.g.dir/sample_name).")         val output: PathPrefix,
  @arg(flag='S', doc="The name of the sample.")                           val sample: String,
  @arg(flag='L', doc="The name of the library.")                          val library: String,
  @arg(flag='P', doc="The platform unit (@RG.PU).  Either one value or one per pair of FASTQs.") val platformUnit: List[String],
  @arg(flag='t', doc="Path to a temporary directory.  Use output if none is given.") val tmp: Option[DirPath] = None,
  @arg(flag='A', doc="The read structure for read one.") val readStructureReadOne: String,
  @arg(flag='B', doc="The read structure for read two.") val readStructureReadTwo: String,
  @arg(flag='r', doc="Path to the reference FASTA.")              val ref: PathToFasta,
  @arg(flag='l', doc="Regions to analyze.")                       val intervals: PathToIntervals,
  @arg(flag='u', doc="The tag containing the raw UMI.")           val umiTag: String = "RX",
  @arg(flag='m', doc="Minimum mapping quality to include reads.") val minMapQ: Int = 10,
  @arg(flag='x', doc="The allowable number of edits between UMIs.") val edits: Int = 1,
  @arg(flag='M', minElements=1, maxElements=3, doc="The minimum number of raw reads per source molecule.")
  val minReads: Seq[Int] = Seq(1),
  @arg(flag='s', doc="Call genotypes per-duplex-strand") val perStrand: Boolean = false,
  @arg(          doc="Keep intermediate files when genotyping.") val keepIntermediates: Boolean = false

) extends Pipeline(outputDirectory=Some(output.getParent)) {

  def build(): Unit = {

    val prepare = new PrepareUnmappedBamPipeline(
      fastq1               = fastq1,
      fastq2               = fastq2,
      output               = output,
      sample               = sample,
      library              = library,
      platformUnit         = platformUnit,
      tmp                  = tmp,
      readStructureReadOne = readStructureReadOne,
      readStructureReadTwo = readStructureReadTwo
    )

    val mapAndGroup = new MapAndGroupRawReads(
      unmappedBam = prepare.unmappedBamFile,
      ref = ref,
      intervals = intervals,
      output = output,
      umiTag = umiTag,
      minMapQ = minMapQ,
      edits = edits
    )

    val genotype = new GenotypeFromGroupedBam(
      input = mapAndGroup.groupedBamFile,
      ref = ref,
      intervals = intervals,
      output = output,
      minReads = minReads,
      perStrand = perStrand,
      keepIntermediates = keepIntermediates,
      tmp = tmp
    )

    root ==> prepare ==> mapAndGroup ==> genotype

  }
}