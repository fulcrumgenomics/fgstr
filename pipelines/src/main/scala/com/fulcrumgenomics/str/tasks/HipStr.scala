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

package com.fulcrumgenomics.str.tasks

import java.nio.file.Path

import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.str.tasks.HipStr.PathToBed
import dagr.core.config.Configuration
import dagr.core.execsystem.{Cores, Memory}
import dagr.core.tasksystem.{FixedResources, ProcessTask, SimpleInJvmTask}
import dagr.tasks.DagrDef.{FilePath, PathToBam, PathToFasta, PathToVcf}
import com.fulcrumgenomics.commons.CommonsDef._
import htsjdk.samtools.util.IntervalList

import scala.collection.mutable.ListBuffer

object HipStr extends Configuration {
  type PathToBed = FilePath
  val HipStrExecutableConfigKey: String = "hipstr.executable"

  def findHipStr: Path = configureExecutableFromBinDirectory(HipStrExecutableConfigKey, "HipSTR")

   /** Creates the regions BED file for HipStr: Required format is tab-delimited columns CHROM START STOP PERIOD NCOPIES */
  class CreateRegionsBed(intervals: PathToIntervals, bed: PathToBed) extends SimpleInJvmTask {
    def run(): Unit = {
      val list  = IntervalList.fromFile(intervals.toFile)
      val lines = list.map { interval =>
        val nameTokens = interval.getName.split(',')
        require(nameTokens.length == 3 || nameTokens.length == 5, s"Require 3 or 5 fields in the name, found ${nameTokens.length} for interval: $interval")
        Seq(
          interval.getContig,
          interval.getStart-1,
          interval.getEnd,
          nameTokens(0),
          nameTokens(1),
          nameTokens(2)
        ).mkString("\t")
      }.toSeq
      Io.writeLines(path=bed, lines=lines)
    }
  }
}

class HipStr(input: PathToBam,
             ref: PathToFasta,
             regions: PathToBed,
             output: PathToVcf,
             stutterIn: Option[FilePath] = None,
             stutterOut: Option[FilePath] = None,
             genotypeLikelihoods: Boolean = true,
             posteriorLikelihoods: Boolean = true,
             minReads: Option[Int] = Some(1),
             useUnpaired: Boolean = false,
             removeDuplicates: Boolean = false,
             haploidChromosomes: Option[IntervalList] = None,
             maxStrLength: Option[Int] = Some(150),
             maxFlankHaplotypes: Option[Int] = Some(10),
             filterFlankingHaplotypes: Boolean = true)
  extends ProcessTask with FixedResources {

  requires(Cores(1), Memory("2g"))

  // Find it when building so we fail earlier
  private val hipStr = HipStr.findHipStr

  override def args: Seq[Any] = {
    val buffer = ListBuffer[Any]()

    buffer.append(hipStr)
    buffer.append("--bams", input)
    buffer.append("--fasta", ref)
    buffer.append("--regions", regions)
    buffer.append("--str-vcf", output)
    stutterIn.foreach(s => buffer.append("--stutter-in", s))
    stutterOut.foreach(s => buffer.append("--stutter-out", s))
    if (genotypeLikelihoods) buffer.append("--output-gls")
    if (posteriorLikelihoods) buffer.append("--output-pls")
    minReads.foreach(m => buffer.append("--min-reads", m))
    if (useUnpaired) buffer.append("--use-unpaired")
    if (!removeDuplicates) buffer.append("--no-rmdup")
    haploidChromosomes.foreach { intervals =>
      if (intervals.nonEmpty) {
        buffer.append("--haploid-chrs")
        buffer.append(intervals.map(_.getContig).toSeq.distinct.mkString(","))
      }
    }
    maxStrLength.foreach(l => buffer.append("--max-str-len", l))
    maxFlankHaplotypes.foreach(l => buffer.append("--max-flank-haps", l))
    if (filterFlankingHaplotypes) buffer.append("--filter-flank-haps")

    buffer
  }
}