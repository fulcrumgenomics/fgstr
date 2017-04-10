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

import dagr.commons.io.Io
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.{Pipeline, ShellCommand}
import dagr.sopt.{arg, clp}
import dagr.tasks.DagrDef._
import dagr.tasks.bwa.Bwa
import dagr.tasks.fgbio._
import dagr.tasks.picard._

@clp(
  description =
    """
      |Pipeline to create duplex-UMI consensus reads and assess them.
      |
      |Currently, maps the reads with BWA, marks duplicates, groups the reads by UMI to assign each read to a source
      |molecule, calls duplex consensus reads, re-maps the duplex consensus reads, filters them, and calls variants.
      |
      |Metrics are then generated on each of the RAW aligned bam, RAW deduped BAM and duplex consensus BAM.
    """,
  group = classOf[Pipelines]
)
class DuplexConsensusPipeline
(@arg(          doc="Path to the unmapped BAM file.")            val unmappedBam: PathToBam,
 @arg(flag="r", doc="Path to the reference FASTA.")              val ref: PathToFasta,
 @arg(flag="l", doc="Regions to analyze.")                       val intervals: PathToIntervals,
 @arg(flag="o", doc="Path prefix for output files.")             val output: PathPrefix,
 @arg(flag="u", doc="The tag containing the raw UMI.")           val umiTag: String = "RX",
 @arg(flag="m", doc="Minimum mapping quality to include reads.") val minMapQ: Int = 10,
 @arg(flag="e", doc="The allowable number of edits between UMIs.") val edits: Int = 1,
 @arg(flag="1", doc="The Phred-scaled error rate for an error prior to the UMIs being integrated.") val errorRatePreUmi: Option[Int] = Some(45),
 @arg(flag="2", doc="The Phred-scaled error rate for an error post the UMIs have been integrated.") val errorRatePostUmi: Option[Int] = Some(30),
 @arg(          doc="The minimum input base quality for bases to go into consensus.")               val minInputBaseQuality: Option[Int] = Some(20),
 @arg(flag="N", doc="Mask (make 'N') consensus bases with quality less than this threshold.")       val minConsensusBaseQuality: Option[Int] = Some(50),
 @arg(flag="M", minElements=1, maxElements=3, doc="The minimum number of reads supporting a consensus base/read.")
 val minReads: Seq[Int] = Seq(6,3,3),
 @arg(flag="E", minElements=1, maxElements=3, doc="The maximum raw-read error rate across the entire consensus read.")
 val maxReadErrorRate: Seq[Double] = Seq(0.05),
 @arg(flag="e", minElements=1, maxElements=3, doc="The maximum error rate for a single consensus base.")
 val maxBaseErrorRate: Seq[Double] = Seq(0.2),
 @arg(flag="n", doc="Maximum fraction of no-calls in the read after filtering.")
 val maxNoCallFraction: Double = 0.2,
 @arg(flag="f", doc="The minimum allele frequency for variant calling") val minimumAf: Double = 0.0001
)
  extends Pipeline(Some(output.getParent)) {

  Io.assertReadable(unmappedBam)
  Io.assertCanWriteFile(output, parentMustExist=false)

  override def build(): Unit = {
    // Files that we're going to create
    val dir = output.getParent
    val pre = output.getFileName
    val mappedRaw   = dir.resolve(pre + ".raw.aligned.bam")
    val dedupedRaw  = dir.resolve(pre + ".raw.deduped.bam")
    val grouped     = dir.resolve(pre + ".grouped.bam")
    val familySizes = dir.resolve(pre + ".family_sizes.txt")
    val dsUnmapped  = dir.resolve(pre + ".consensus.unmapped.bam")
    val dsMapped    = dir.resolve(pre + ".consensus.mapped.bam")

    // Copy the intervals to the directory and make a BED file
    val regionsIl     = dir.resolve("calling_regions.interval_list")
    val copyIntervals = new ShellCommand("cp", this.intervals.toString, regionsIl.toString)

    // Make the grouped and consensus BAMs
    val bwaRaw    = Bwa.bwaMemStreamed(unmappedBam=unmappedBam, mappedBam=mappedRaw, ref=ref, samToFastqCores=0, bwaMemMemory="6G")
    val group     = new GroupReadsByUmi(in=mappedRaw, out=grouped, familySizeOut=Some(familySizes),
      rawTag=Some(umiTag), minMapQ=Some(minMapQ), strategy=AssignmentStrategy.Paired, edits=Some(edits))

    val callConsensus = new CallDuplexConsensusReads(
      in=grouped,
      out=dsUnmapped,
      readNamePrefix=None,
      errorRatePreUmi=errorRatePreUmi,
      errorRatePostUmi=errorRatePostUmi,
      minInputBaseQuality=minInputBaseQuality
    )

    val bwaDs = Bwa.bwaMemStreamed(unmappedBam=dsUnmapped, mappedBam=dsMapped, ref=ref, samToFastqCores=0, bwaMemMemory="6G", mergeBamAlignmentMem="6G")

    val postPipeline = new PostDuplexConsensusPipeline(input=dir, ref=ref, intervals=regionsIl, output=dir,
      minConsensusBaseQuality=minConsensusBaseQuality, minReads=minReads, maxBaseErrorRate=maxBaseErrorRate,
      maxNoCallFraction=maxNoCallFraction, minimumAf=minimumAf)

    root ==> bwaRaw ==> group ==> callConsensus ==> bwaDs ==> (postPipeline :: new DeleteBam(mappedRaw, dsUnmapped))
    root ==> copyIntervals ==> postPipeline
  }
}
