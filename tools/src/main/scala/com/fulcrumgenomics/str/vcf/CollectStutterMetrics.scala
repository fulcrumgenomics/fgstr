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

package com.fulcrumgenomics.str.vcf

import com.fulcrumgenomics.FgBioDef.{PathPrefix, PathToIntervals, PathToVcf}
import com.fulcrumgenomics.cmdline.ClpGroups
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.str.cmdline.FgStrTool
import htsjdk.variant.vcf.VCFFileReader
import com.fulcrumgenomics.commons.CommonsDef.SafelyClosable
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.str.vcf.StrInterval.StrAllele
import com.fulcrumgenomics.util.Rscript
import htsjdk.variant.variantcontext.VariantContext

import scala.collection.mutable
import scala.util.Failure

@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Collects the histogram of stutter events.
    |
    |## Inputs
    |
    |The VCF from `StrGenotypeDuplexMolecules` should be given as input.  Optionally, one or more VCFs produced by
    |`HipSTR` can also be given.  The VCFs from `HipSTR` should have disjoint variant calls (ex. one per STR).
    |
    |An interval list specifying the set of regions over which to call STRs should be given. The name field should
    |contain a comma list of values as follows:
    |  1. the repeat unit length (ex. `3` for the tri-nucleotide repeat `TCATCATCATCA`).
    |  2. the number of repeat units (ex. `4` for the tri-nucleotide repeat `TCATCATCATCA`).
    |  3. the name of the STR (ex. D1S1656)
    |Optionally, additional columns can be given for one or more expected truth alleles.  For example, a known haploid
    |call should have one extra column, a known diploid call should have two extra columns, and so on.
    |
    |## Outputs
    |
    |  * `<output>.stutter.tab> - attempts to estimate the frequency of stutter events.  A stutter occurs when we have
    |    one or more bases missing or present relative to the ground truth.  The stutter distance is the minimum # of
    |    bases different between the given call and any ground truth allele.  If no truth alleles are given, the
    |    called alleles are used instead.  The output will include a row per stutter distance (ex. -1, 0, or 1) and a
    |    column with counts per STR site, as well as a final column with counts aggregated across all sites.
    |  * `<output.stutter.raw.tab> - the same as the previous but computed using the raw read counts.
    |  * `<output>.stutter.pdf` - a plot per STR with the histogram of counts across stutter distances.
  """)
class CollectStutterMetrics
(
  @arg(flag='i', doc="Input VCF of variant calls from `StrGenotypeDuplexMolecules`.") val input: PathToVcf,
  @arg(flag='l', doc="Interval list with the STR regions and optionally known calls.") val intervals: PathToIntervals,
  @arg(flag='o', doc="Prefix for all output files.") val output: PathPrefix,
  @arg(flag='r', doc="Input VCF(s) of variant calls from `HipSTR`") val hipStrCalls: Seq[PathToVcf],
  @arg(flag='m', doc="Minimum % of 0-stutter observations to include in the aggregate counts") val minZeroStutter: Double = 0.5,
  private val skipPlots: Boolean = false // for not plotting in tests
) extends FgStrTool with LazyLogging {

  private val ScriptPath = "com/fulcrumgenomics/str/vcf/CollectStutterMetrics.R"

  Io.assertReadable(input)
  Io.assertReadable(intervals)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val vcfIn          = new VCFFileReader(input.toFile, false)
    val vcfHipStr      = hipStrCalls.map(in => new VCFFileReader(in.toFile, false))
    val strIntervals   = StrInterval.loadIntervals(this.intervals).toSeq
    val counters       = new StutterCounter()
    val hipStrCounters = new StutterCounter()

    // Get the counts per STR
    strIntervals.foreach { str =>
      // get the "ref" calls from the genotyped bam if not given in the interval list
      val refCalls: Seq[Float] = vcfIn.query(str.chrom, str.start, str.end)
        .filter(ctx => ctx.getNoCallCount != ctx.getNSamples)
        .toSeq match {
        case Seq(ctx) =>
          str.truthCalls match {
            case Seq()      => ctx.getGenotype(0).getAttributeAsString("STR_GT", "").split(",").map(_.toFloat)
            case truthCalls => truthCalls
          }
        case calls =>
          Seq.empty
      }

      if (refCalls.isEmpty) {
        counters.add(str.name, new NumericCounter[Int]())
        hipStrCounters.add(str.name, new NumericCounter[Int]())
      }
      else {
        // genotyped bam
        counters.add(str.name, process(str, vcfIn, refCalls))

        // raw reads
        vcfHipStr.foreach { vcf =>
          hipStrCounters.add(str.name, process(str, vcf, refCalls))
        }
      }
    }

    // genotyped bam
    val stutterOut = PathUtil.pathTo(output + ".stutter.tab")
    counters.write(PathUtil.pathTo(output + ".stutter.tab"), strIntervals, minZeroStutter)

    // raw reads
    val rawStutterOut = PathUtil.pathTo(output + ".stutter.raw.tab")
    hipStrCounters.write(rawStutterOut, strIntervals, minZeroStutter)

    vcfIn.safelyClose()

    if (!skipPlots) {
      logger.info("Plotting")
      val pdf = PathUtil.pathTo(output + ".stutter.pdf")
      Rscript.execIfAvailable(ScriptPath, stutterOut.toString, rawStutterOut.toString, pdf.toString) match {
        case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
        case _ => Unit
      }
    }
    else {
      logger.warning("Skipping plots")
    }
  }

  // find the one closes to zero
  private def minAbs(x: Double*): Double = {
    val min = x.map(Math.abs).min
    x.find { y => Math.abs(y) == min }.get
  }

  /** Gets the calls for this STR and returns a [[NumericCounter]] that counts the # of observations for
    * each stutter event. */
  private def process(str: StrInterval, vcfIn: VCFFileReader, refCalls: Seq[Float]): NumericCounter[Int] = {
    vcfIn.query(str.chrom, str.start, str.end)
      .filter(ctx => ctx.getNoCallCount != ctx.getNSamples)
      .toSeq match {
        case Seq(ctx) => if (ctx.getNSamples == 1) collectSingle(str, ctx, refCalls) else collectRaw(str, ctx, refCalls)
        case calls => new NumericCounter[Int]()
      }
  }

  /** Given an STR and set of calls, returns a [[NumericCounter]] that counts the # of observations for
    * each stutter event. */
  private def collect(str: StrInterval, refCalls: Seq[Float], calls: Seq[StrAllele]): NumericCounter[Int] = {
    val counter = new NumericCounter[Int]()
    calls.filter(_.count > 0)
      .foreach { call =>
        val stutter = {
          val abs = minAbs(refCalls.map(refCall => call.repeatLength - refCall):_*)
          if (abs < 0) Math.ceil(abs) else Math.floor(abs)
        }.toInt
        counter.count(stutter, call.count)
      }
    counter
  }

  /** Assumes a variant context produced by `StrGenotypeDuplexMolecules` and returns a [[NumericCounter]] that counts
    * the # of observations for each stutter event. */
  private def collectSingle(str: StrInterval, ctx: VariantContext, refCalls: Seq[Float]): NumericCounter[Int] = {
    // Get all calls seen, not just the one genotyped
    val counts = {
      val refCount  = ctx.getAttributeAsInt("REFAC", 0)
      val altCounts = ctx.getAttributeAsIntList("AC", 0).map(_.toInt).toSeq
      refCount +: altCounts
    }
    val calls = StrAllele.toCalls(str, ctx, counts)

    collect(str, refCalls, calls)
  }

  /** Assumes a variant context produced by `HipSTR` and returns a [[NumericCounter]] that counts
    * the # of observations for each stutter event. */
  private def collectRaw(str: StrInterval, ctx: VariantContext, refCalls: Seq[Float]): NumericCounter[Int] = {
    val counter = new NumericCounter[Int]()

    ctx.getGenotypes.filterNot(_.isNoCall).foreach { genotype =>
      require(genotype.getAlleles.length == 1, s"Expected a single allele, found '${genotype.getAlleles.length}'")

      val calls: Seq[StrAllele] = genotype.getAttributeAsString("ALLREADS", "")
        .split(';')
        .filter(_.nonEmpty)
        .map { keyVal =>
          keyVal.split('|').toSeq match {
            case Seq(bpDiff, count) =>
              // NB: bpDiff is the # of bases different relative to the reference!
              val alleleLength = bpDiff.toInt + (str.refLength * str.unitLength)
              StrAllele(allele=null, alleleLength=alleleLength, count=count.toInt, unitLength=str.unitLength)
            case _ => throw new IllegalArgumentException(s"Could not parse ALLREADS entry: $keyVal")
          }
       }

      collect(str, refCalls, calls).iterator.foreach { case (repeatLength, count) =>
         counter.count(repeatLength, count)
      }
    }

    counter
  }
}

private object StutterCounter {
  val AllStrsKey: String = "all_strs"
}

private class StutterCounter(stutters: Seq[Int] = Seq.range(-5,6)) {
  private val counters = mutable.HashMap[String,NumericCounter[Int]]()

  def add(name: String, counter: NumericCounter[Int]): Unit = {
    if (counters.contains(name)) {
      val strCounter = counters(name)
      counter.iterator.foreach { case (stutter, count) => strCounter.count(stutter, count) }
    }
    else {
      counters(name) = counter
    }
  }

  def write(path: FilePath, strIntervals: Seq[StrInterval], minZeroStutter: Double): Unit = {
    // Aggregate the counts across all STRs
    val allStrCounter = new NumericCounter[Int]()
    counters.foreach { case (_, counter) =>
      // Do not include sites where the fraction of reads/calls for zero stutters is too low.  This may be caused by
      // having the wrong genotype.
      val totalCounts = counter.map(_._2).sum
      val zeroStutterFreq = counter.countOf(0) / totalCounts.toDouble
      if (minZeroStutter <= zeroStutterFreq) {
        counter.foreach { case (stutter, count) => allStrCounter.count(stutter, count) }
      }
    }
    counters(StutterCounter.AllStrsKey) = allStrCounter

    // Print the metrics output
    val names    = strIntervals.map(_.name) ++ Seq("all_strs")
    val header   = (Seq("STR") ++ stutters).mkString("\t")
    val lines = names.map { name =>
      val counter = counters(name)
      val values = stutters.map { stutter => counter.countOf(stutter) }
      (Seq(name) ++ values).mkString("\t")
    }
    Io.writeLines(path, Seq(header) ++ lines)
  }
}