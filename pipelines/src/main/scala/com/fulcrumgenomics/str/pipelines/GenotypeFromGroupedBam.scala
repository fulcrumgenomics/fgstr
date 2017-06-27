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

import java.nio.file.Files

import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.str.tasks.HipStr.CreateRegionsBed
import com.fulcrumgenomics.str.tasks.{HipStr, ReadGroupPerDuplexMolecularId, StrGenotypeDuplexMolecules}
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.{Pipeline, SimpleInJvmTask}
import dagr.tasks.DagrDef.{FilePath, PathPrefix, PathToBam, PathToFasta, PathToIntervals}
import dagr.tasks.ScatterGather.{Partitioner, Scatter}
import dagr.tasks.misc.{DeleteFiles, IndexVcfGz}
import dagr.tasks.picard.GatherVcfs
import htsjdk.samtools.util.IntervalList

@clp(
  description =
    """
      |Pipeline to call STR genotypes from a grouped BAM.
      |
      |Assigns reads from the same source molecule to their own read group, filters the reads from the same source
      |duplex molecule to require a minimum number of reads (overall, AB strand, and BA strand), runs HipSTR to call
      |genotypes jointly across all source duplex molecules (treats each source molecule as its own sample), and
      |combines the per-source-duplex-molecule genotypes into a final genotype.
      |
      |The grouped BAM should be created by fgbio's `GroupReadsByUmi` with the `--strategy=paired` option or by the
      |`MapAndGroupRawReads` pipeline.
      |
      |An interval list specifying the set of regions over which to call STRs should be given. The name field should
      |contain a comma list of values as follows:
      |  1. the repeat unit length (ex. `3` for the tri-nucleotide repeat `TCATCATCATCA`).
      |  2. the number of repeat units (ex. `4` for the tri-nucleotide repeat `TCATCATCATCA`).
      |  3. the name of the STR (ex. D1S1656)
      |  4. optionally, the expected (known or truth) number of repeat units for allele #1
      |  5. optionally, the expected (known or truth) number of repeat units for allele #2
      |An example name field with the optionals is `4,17,D1S1656,9,10` and without is `4,17,D1S1656`.
      |
      |The option `--min-reads` make take take 1-3 values. For example:
      |
      |```
      |FilterRawReadsFromDuplexMolecule ... --min-reads 10 5 3
      |```
      |
      |In each case if fewer than three values are supplied, the last value is repeated (i.e. `80 40` -> `80 40 40`
      |and `10` -> `10 10 10`.  The first value applies to all the reads from the same source molecule, the second value
      |to one single-strand consensus, and the last value to the other single-strand consensus. It is required that if
      |values two and three differ, the _more stringent value comes earlier_.
      |
      |The `--per-strand` option will call group reads by duplex strand (i.e. A and B strand respectively).  In this case,
      |HipSTR will be called twice: once to learn the stutter model by calling genotypes treating reads from the same
      |duplex source molecule, and a second time using the stutter model to call genotypes treating reads from the same
      |duplex source molecule and strand.  The genotype calls for a given duplex source molecule (A and B strand respectively)
      |are compared for perfect agreement to create a per duplex source molecule genotype call, which then is used as
      |describe above.
    """,
  group = classOf[Pipelines]
)
class GenotypeFromGroupedBam
(@arg(flag='i', doc="Input SAM or BAM file (from GroupReadsByUmi).") val input: PathToBam,
 @arg(flag='r', doc="Path to the reference FASTA.")                  val ref: PathToFasta,
 @arg(flag='l', doc="STR regions to analyze.")                       val intervals: PathToIntervals,
 @arg(flag='o', doc="Path prefix for output files.")                 val output: PathPrefix,
 @arg(flag='M', minElements=1, maxElements=3, doc="The minimum number of raw reads per source molecule.")
 val minReads: Seq[Int] = Seq(1),
 @arg(flag='s', doc="Call genotypes per-duplex-strand") val perStrand: Boolean = false,
 @arg(flag='t', doc="Temporary directory in which to store intermediate results.") val tmp: Option[DirPath] = None
) extends Pipeline(Some(output.getParent)) {


  override def build(): Unit = {
    Io.assertReadable(Seq(input, ref, intervals))
    Io.assertCanWriteFile(output)
    tmp.foreach(Io.assertListable)

     def toStrName(intervalList: IntervalList): String = {
      require(intervalList.length == 1, s"Expected one STR interval, found ${intervalList.length}")
      val interval = intervalList.iterator.next
      val tokens   = interval.getName.split(',')
      require(3 <= tokens.length & tokens.length <= 5, s"Expected name for the interval to have 3-5 comma-separated fields: $interval")
      tokens(2) // the name of the STR
    }

    val dir = tmp match {
      case Some(p) => Files.createTempDirectory(p, "genotype_from_grouped_bam.")
      case None    => Files.createTempDirectory("genotype_from_grouped_bam.")
    }
    val scatterer = new SplitStrIntervals(intervals=this.intervals)
    val scatter   = Scatter(scatterer)
    val genotypes = scatter.map { intervalList =>
      val strName = toStrName(intervalList)
      new GenotypeStr(
        input        = input,
        ref          = ref,
        intervalList = intervalList,
        minReads     = minReads,
        perStrand    = perStrand,
        output       = dir.resolve(strName)
      )
    }
    val gather = genotypes.gather { tasks: Seq[GenotypeStr] =>
      new GatherVcfs(
        in  = tasks.map { task => PathUtil.pathTo(task.output + ".vcf.gz") },
        out = PathUtil.pathTo(this.output + ".vcf.gz")
      )
    }
    root ==> scatter
    gather ==> new DeleteFiles(dir)
  }
}

/** A pipeline that calls a genotype for a single STR. This pipeline does all the heavy lifting! */
private class GenotypeStr
( val input: PathToBam,
  val ref: PathToFasta,
  val intervalList: IntervalList,
  val minReads: Seq[Int] = Seq(1),
  val perStrand: Boolean = false,
  val output: DirPath,
  val suffix: Option[String] = None
) extends Pipeline(suffix=Some("." + output.getFileName)) {

  def f(ext: String): FilePath = PathUtil.pathTo(output + ext)

  require(intervalList.length == 1, s"Expected one STR interval, found ${intervalList.length}")

  def build(): Unit = {
    val bed       = f(".regions.bed")
    val midBam    = f(".mid.bam")
    val hipStrVcf = f(".hipstr.vcf.gz")
    val stutter   = f(".stutter.txt")
    val intervals = f(".interval_list")

    val toIntervals    = new SimpleInJvmTask {
      name = "WriteStrInterval"
      override def run(): Unit = intervalList.write(intervals.toFile)
    }
    val toBed          = new CreateRegionsBed(intervals=intervals, bed=bed)

    val toHipstrVcf = if (perStrand) {
      val strandBam = f("mid.per_strand.bam")
      // learn the stutter model from calling reads from the same duplex source molecule
      val toDuplexRgBam  = new ReadGroupPerDuplexMolecularId(input=input, output=midBam, ref=Some(ref), intervals=Some(intervals), minReads=minReads)
      val toStutterModel = new HipStr(input=midBam, ref=ref, regions=bed, output=hipStrVcf, stutterOut=Some(stutter), haploidChromosomes=Some(intervalList))

      // apply the stutter model when calling reads from the same strand of the duplex source molecule
      val toStrandRgBam  = new ReadGroupPerDuplexMolecularId(input=input, output=strandBam, ref=Some(ref), intervals=Some(intervals), minReads=minReads, perStrand=true)
      val toHipstrVcf    = new HipStr(input=strandBam, ref=ref, regions=bed, output=hipStrVcf, stutterIn=Some(stutter), haploidChromosomes=Some(intervalList))

      toIntervals ==> (((toBed :: toDuplexRgBam) ==> toStutterModel) :: toStrandRgBam) ==> toHipstrVcf
    }
    else {
      val toRgBam     = new ReadGroupPerDuplexMolecularId(input=input, output=midBam, ref=Some(ref), intervals=Some(intervals), minReads=minReads)
      val toHipstrVcf = new HipStr(input=midBam, ref=ref, regions=bed, output=hipStrVcf, haploidChromosomes=Some(intervalList))
      toIntervals ==> (toBed :: toRgBam) ==> toHipstrVcf
    }

    // merge the calls
    val indexVcf      = new IndexVcfGz(hipStrVcf)
    val toFinalVcf    = new StrGenotypeDuplexMolecules(input=hipStrVcf, output=output, ref=ref, intervals=intervals, perStrand=perStrand)
    root ==> toHipstrVcf ==> indexVcf ==> toFinalVcf
  }
}

/** Splits the given interval list into one file per interval. */
private class SplitStrIntervals(intervals: PathToIntervals)
  extends SimpleInJvmTask with Partitioner[IntervalList] {

  private var intervalLists : Option[Seq[IntervalList]] = None

  override def partitions: Option[Seq[IntervalList]] = this.intervalLists

  def run(): Unit = {
    val intervals = IntervalList.fromFile(this.intervals.toFile)
    val lists     = intervals.zipWithIndex.map { case (interval, idx) =>
      val intervalList = new IntervalList(intervals.getHeader)
      intervalList.add(interval)
      intervalList
    }
    this.intervalLists = Some(lists.toSeq)
  }
}
