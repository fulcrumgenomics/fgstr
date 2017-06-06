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

import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.commons.io.Io
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.{Pipeline, ShellCommand}
import dagr.tasks.DagrDef._
import dagr.tasks.bwa.Bwa
import dagr.tasks.fgbio._
import dagr.tasks.picard._

@clp(
  description =
    """
      |Pipeline to map and group raw reads.
      |
      |Currently, maps the raw reads with BWA, groups the reads by UMI to assign each read to a source
      |molecule.
      |
      |Metrics are then generated on the raw aligned bam.
    """,
  group = classOf[Pipelines]
)
class MapAndGroupRawReads
(@arg(          doc="Path to the unmapped BAM file.")            val unmappedBam: PathToBam,
 @arg(flag='r', doc="Path to the reference FASTA.")              val ref: PathToFasta,
 @arg(flag='l', doc="Regions to analyze.")                       val intervals: PathToIntervals,
 @arg(flag='o', doc="Path prefix for output files.")             val output: PathPrefix,
 @arg(flag='u', doc="The tag containing the raw UMI.")           val umiTag: String = "RX",
 @arg(flag='m', doc="Minimum mapping quality to include reads.") val minMapQ: Int = 10,
 @arg(flag='x', doc="The allowable number of edits between UMIs.") val edits: Int = 1
)
  extends Pipeline(Some(output.getParent)) {

  override def build(): Unit = {

    Io.assertReadable(unmappedBam)
    Io.assertCanWriteFile(output, parentMustExist=false)

    // Files that we're going to create
    val dir = output.getParent
    val pre = output.getFileName
    val mappedRaw   = dir.resolve(pre + ".raw.aligned.bam")
    val grouped     = dir.resolve(pre + ".grouped.bam")
    val familySizes = dir.resolve(pre + ".family_sizes.txt")

    // Make the grouped and consensus BAMs
    val bwaRaw    = Bwa.bwaMemStreamed(unmappedBam=unmappedBam, mappedBam=mappedRaw, ref=ref, samToFastqCores=0, bwaMemMemory="6G")
    val group     = new GroupReadsByUmi(in=mappedRaw, out=grouped, familySizeOut=Some(familySizes),
      rawTag=Some(umiTag), minMapQ=Some(minMapQ), strategy=AssignmentStrategy.Paired, edits=Some(edits))
    val dsMetrics         = new CollectDuplexSeqMetrics(input=grouped, output=dir.resolve(pre), intervals=Some(intervals))
    val dsMetricsAllReads = new CollectDuplexSeqMetrics(input=grouped, output=dir.resolve(pre + ".all_reads"), intervals=None)

    root ==> bwaRaw ==> group ==> (dsMetrics :: dsMetricsAllReads :: new DeleteBam(mappedRaw))
  }
}
