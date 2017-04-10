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

import dagr.commons.io.{Io, PathUtil}
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.Pipeline
import dagr.sopt.{arg, clp}
import dagr.tasks.DagrDef.{PathPrefix, PathToBam, PathToFasta, PathToIntervals, PathToVcf}
import dagr.tasks.fgbio.{ClipOverlappingReads, FilterConsensusReads}
import dagr.tasks.picard.IntervalListToBed
import dagr.tasks.vc.VarDictJavaEndToEnd
import htsjdk.variant.vcf.VCFFileReader


@clp(
  description =
    """
      |Pipeline to filter duplex consensus reads and compute metrics.
      |
      |Requires as input the output path prefix given to DuplexConsensusPipeline.
    """,
  group = classOf[Pipelines]
)
class PostDuplexConsensusPipeline
(
  @arg(flag="i", doc="Path to the input consensus duplex BAM.") val input: PathToBam,
  @arg(flag="r", doc="Path to the reference FASTA.")            val ref: PathToFasta,
  @arg(flag="l", doc="Regions to analyze.")                     val intervals: PathToIntervals,
  @arg(flag="o", doc="Path prefix for output files.")           val output: PathPrefix,
  @arg(flag="N", doc="Mask (make 'N') consensus bases with quality less than this threshold.")
  val minConsensusBaseQuality: Option[Int] = Some(50),
  @arg(flag="M", minElements=1, maxElements=3, doc="The minimum number of reads supporting a consensus base/read.")
  val minReads: Seq[Int] = Seq(6,3,3),
  @arg(flag="E", minElements=1, maxElements=3, doc="The maximum raw-read error rate across the entire consensus read.")
  val maxReadErrorRate: Seq[Double] = Seq(0.05),
  @arg(flag="e", minElements=1, maxElements=3, doc="The maximum error rate for a single consensus base.")
  val maxBaseErrorRate: Seq[Double] = Seq(0.1),
  @arg(flag="n", doc="Maximum fraction of no-calls in the read after filtering.")
  val maxNoCallFraction: Double = 0.2,
  @arg(flag="f", doc="The minimum allele frequency for variant calling") val minimumAf: Double = 0.0001
)
  extends Pipeline(Some(output.getParent)) {

  def build(): Unit = {

    Io.assertReadable(input)
    Io.assertCanWriteFile(output, parentMustExist=false)

    // Files that we're going to create
    val outDir      = output.getParent
    val inDir       = input.getParent
    val dsFiltered  = outDir.resolve(output.getFileName + ".consensus.filtered.bam")

    // Make a BED file for the calling region.
    val regionsBed  = inDir.resolve("calling_regions.bed")
    val makeBed     = new IntervalListToBed(intervals=intervals, bed=regionsBed)

    // Filter the duplex consensus reads
    val filterDs = FilterConsensusReads(in=input, out=dsFiltered, ref=ref, minReads=Seq(6,3,3), maxReadErrorRate=Seq(0.05),
      maxBaseErrorRate=Seq(0.1), minQuality=minConsensusBaseQuality.getOrElse(0), maxNoCallFraction=0.05).requires(cores=1, memory="8G")

    // Make variant calls and assess the results
    val clippedBam      = PathUtil.replaceExtension(input, ".clipped.bam")
    val clip            = new ClipOverlappingReads(in=input, out=clippedBam, ref=ref).requires(cores=1, memory="8G")
    val metricsPipeline = new CollectMetricsPipeline(input=clippedBam, ref=ref, intervals=intervals, output=output)

    // Do some variant calling
    val vcf = PathUtil.replaceExtension(input, ".vcf")
    val vardict =  new VarDictJavaEndToEnd(
      tumorBam       = clippedBam,
      tumorName      = Some("tumor"),
      bed            = regionsBed,
      ref            = ref,
      out            = vcf,
      minimumAf      = minimumAf,
      minimumQuality = Some(30),
      pileupMode     = true
    )

    root ==> (filterDs :: makeBed)
    filterDs ==> clip ==> metricsPipeline
    (clip :: makeBed) ==> vardict
  }
}
