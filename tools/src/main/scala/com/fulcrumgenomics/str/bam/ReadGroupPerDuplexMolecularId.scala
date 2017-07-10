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
package com.fulcrumgenomics.str.bam

import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.ClpGroups
import com.fulcrumgenomics.commons.CommonsDef.IteratorToJavaCollectionsAdapter
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.str.FgStrDef
import com.fulcrumgenomics.str.FgStrDef.DuplexFilters
import com.fulcrumgenomics.str.cmdline.FgStrTool
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector}
import htsjdk.samtools.{SAMFileHeader, SAMReadGroupRecord}

import scala.collection.mutable

@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Creates a BAM file with a read group per duplex molecule ID.
    |
    |The output will be coordinate sorted.
    |
    |The raw reads should have been grouped using fgbio's GroupReadsByUmi.
    |
    |The duplex molecular ID for each read should be stored in the given tag (ex. MI).  A new read group per molecular
    |ID will be written to the output BAM file.  The trailing `/A` or `/B` will be removed for duplex reads.
    |
    |The input file should have a single read group, with all fields from this read group copied to each new read group,
    |besides the ID.  The read group ID for reads from a source molecule will be the original read group ID with the
    |source molecule's ID appended (with a dash separator).  The sample name and library will have the molecular ID
    |appended (also with a dash separator).
    |
    |An interval list can be given to filter output reads to a set of intervals.
    |
    |The option `--min-reads` make take take 1-3 values. For example:
    |
    |```
    |ReadGroupPerDuplexMolecularId ... --min-reads 10 5 3
    |```
    |
    |In each case if fewer than three values are supplied, the last value is repeated (i.e. `80 40` -> `80 40 40`
    |and `10` -> `10 10 10`.  The first value applies to all the reads from the same source molecule, the second value
    |to one single-strand consensus, and the last value to the other single-strand consensus. It is required that if
    |values two and three differ, the _more stringent value comes earlier_.
    |
    |This tool will also output the counts per-source molecule of the reads with "/A" and "/B" suffixes on their umi
    |tags to a file with extension ".umi_counts.txt".
  """)
class ReadGroupPerDuplexMolecularId
(
  @arg(flag='i', doc="Input SAM or BAM file (from GroupReadsByUmi).") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM.") val output: PathToBam = Io.StdOut,
  @arg(flag='r', doc="Reference sequence fasta file.") val ref: Option[PathToFasta] = None,
  @arg(flag='l', doc="Optional list of intervals to restrict reads.") val intervals: Option[PathToIntervals] = None,
  @arg(flag='M', minElements=1, maxElements=3, doc="The minimum number of raw reads per source molecule.")
  val minReads: Seq[Int] = Seq(1),
  @arg(flag='T', doc="The output tag from UMI grouping.") val assignTag: String = ConsensusTags.MolecularId,
  @arg(flag='s', doc="Create a read group per-strand of a duplex molecule") val perStrand: Boolean = false,
  @arg(flag='S', doc="The sort order of the output") val samOrder: Option[SamOrder] = None,
  @arg(          doc="Require that reads that span the given intervals (i.e. do not start/stop within)") val span: Boolean = false
) extends FgStrTool with LazyLogging {

  Io.assertReadable(input)
  intervals.foreach(Io.assertReadable)
  Io.assertCanWriteFile(output)

  private val toMolecularId: SamRecord => String = rec => if (this.perStrand) rec[String](assignTag) else FgStrDef.toMolecularId(rec, assignTag).toString
  private val filters = DuplexFilters(minReads)
  private var overlappingReads: Long = 0
  private var numNonSpaningReads: Long = 0

  validate(input.toAbsolutePath != Io.StdIn, "Reading from standard input is not supported")

  override def execute(): Unit = {
    val mids         = new mutable.HashSet[String]()
    val midsByStrand = new SimpleCounter[String]()
    val tmpOut    = PathUtil.replaceExtension(output, ".tmp.bam")
    tmpOut.toFile.deleteOnExit()

    // Pass 1: find all the molecular identifiers to use
    {
      logger.info("Extracting molecular ids")
      val progress = ProgressLogger(this.logger, unit = 5e5.toInt, verb = "Extracted molecule ids")
      val in = SamSource(input, ref = ref)
      val iterator = toIterator(in)
      val tmpWriter = SamWriter(tmpOut, in.header, ref = ref, compression=1)
      iterator.foreach { rec =>
        mids += toMolecularId(rec)
        midsByStrand.count(rec[String](assignTag))
        tmpWriter += rec
        progress.record()
        progress.record(rec)
        if (progress.getCount % 5e5.toInt == 0) {
          logger.info("# of molecules so far: " + mids.iterator.length)
        }
      }
      in.safelyClose()
      tmpWriter.close()
      logger.info(s"Extracted ${mids.iterator.length} molecular ids across ${progress.getCount} raw reads")
      if (span) logger.info(f"Skipped $numNonSpaningReads out of $overlappingReads (${numNonSpaningReads/overlappingReads.toDouble * 100}%.2f%%) overlapping reads due to not fully spanning the STR")
    }

    // Pass 2: update the header and records
    {
      logger.info("Writing to output")
      val in        = SamSource(tmpOut, ref=ref)
      val iterator  = new DuplexMoleculeFilteringIterator(toIterator(in), filters, assignTag)
      val header    = toOutputHeader(in, mids.toSeq.distinct, midsByStrand, samOrder)
      val out       = SamWriter(output, header, ref=ref, sort=samOrder)
      val metricOut = PathUtil.replaceExtension(output, ".umi_counts.txt")
      val passFilter = new mutable.HashSet[String]()

      if (mids.nonEmpty) {
        iterator.foreach { rec =>
          val mid   = toMolecularId(rec)
          rec("RG") = s"mid-$mid"
          passFilter += mid
          out += rec
        }
      }
      in.safelyClose()
      out.close()

      val metrics = mids.toSeq.distinct.map { mid =>
        UmiCounts(
          mid         = mid,
          counts_a    = midsByStrand.countOf(mid + "/A"),
          counts_b    = midsByStrand.countOf(mid + "/B"),
          pass_filter = passFilter.contains(mid)
        )
      }
      Metric.write[UmiCounts](metricOut, metrics)
    }
  }

  /** Create the output header. */
  private def toOutputHeader(in: SamSource, mids: Seq[String], midsByStrand: SimpleCounter[String], sortOrder: Option[SamOrder]): SAMFileHeader = {
    val header = in.header.clone()

    val groupsOfReadGroups: Seq[Seq[SAMReadGroupRecord]] = header.getReadGroups
      .toSeq
      .groupBy { rg => (rg.getSample, rg.getLibrary) }
      .values
      .toSeq

    val sourceReadGroups = groupsOfReadGroups match {
      case Seq()    => throw new IllegalArgumentException(s"Found no read groups in ${this.input}")
      case Seq(rgs) => rgs
      case _        => throw new IllegalArgumentException(s"Expected all read groups to have the same sample/library")
    }

    def toReadGroupValue(f: SAMReadGroupRecord => String, suffix: String = ""): String = {
      val values = sourceReadGroups.map(f).filter(_ != null)
      if (values.isEmpty) null
      else values.distinct.mkString(",") + suffix
    }

    def toReadGroup(mid: String, countA: Long, countB: Long): SAMReadGroupRecord = {
      val readGroup = new SAMReadGroupRecord(s"mid-$mid")
      readGroup.setAttribute("cA", countA.toString) // Custom count
      readGroup.setAttribute("cB", countA.toString) // Custom count
      readGroup.setSample(toReadGroupValue(_.getSample, s"-$mid"))
      readGroup.setLibrary(toReadGroupValue(_.getLibrary, s"-$mid"))
      readGroup.setPlatform(toReadGroupValue(_.getPlatform))
      readGroup.setPlatformModel(toReadGroupValue(_.getPlatformModel))
      readGroup.setPlatformUnit(toReadGroupValue(_.getPlatformUnit))
      readGroup.setDescription(toReadGroupValue(_.getDescription))
      readGroup.setRunDate(if (sourceReadGroups.map(_.getRunDate).distinct.length == 1) sourceReadGroups.head.getRunDate else null)
      readGroup.setSequencingCenter(toReadGroupValue(_.getSequencingCenter))
      readGroup.setPredictedMedianInsertSize(if (sourceReadGroups.map(_.getPredictedMedianInsertSize).distinct.length == 1) sourceReadGroups.head.getPredictedMedianInsertSize else null)
      readGroup
    }

    // create the read groups
    val readGroups = if (mids.nonEmpty) {
      mids.map { mid =>
        if (perStrand) {
          val count = midsByStrand.countOf(mid)
          require(count > 0, s"Bug: counts were zero for mid $mid")
          if (mid.endsWith("/A")) toReadGroup(mid, count, 0) else toReadGroup(mid, 0, count)
        }
        else {
          val countA = midsByStrand.countOf(mid + "/A")
          val countB = midsByStrand.countOf(mid + "/B")
          require(countA + countB > 0, s"Bug: counts were zero for mid $mid")
          toReadGroup(mid, countA, countB)
        }
      }
    }
    else {
      // We may not have ANY mids if there were no reads, or reads did not overlap the intervals,
      // so just create a "dummy" read group.
      Seq(toReadGroup(if (perStrand) "0/A" else "0", 0, 0))
    }

    header.setReadGroups(readGroups.toIterator.toJavaList)
    sortOrder.foreach(_.applyTo(header))
    header
  }

  /** Build the iterator we'll use based on whether or not we're restricting to a set of intervals */
  private def toIterator(in: SamSource): Iterator[SamRecord] = {
    intervals match {
      case None       => in.view.toIterator
      case Some(path) =>
        val ilist = IntervalList.fromFile(path.toFile).uniqued(false)
        val detector = new OverlapDetector[Interval](0, 0)
        detector.addAll(ilist.getIntervals, ilist.getIntervals)
        in.iterator.filter { rec: SamRecord =>
          if (!rec.paired || rec.unmapped || rec.mateUnmapped) false
          else {
            val (start, end) = if (rec.refIndex == rec.mateRefIndex) Bams.insertCoordinates(rec) else (rec.start, rec.end)
            if (span) {
              val overlaps = detector.getOverlaps(new Interval(rec.refName, start, end))
              if (overlaps.nonEmpty) overlappingReads += 1
              // require that at least one end of a pair fully overlaps the str interval
              overlaps.exists { strInterval =>
                val mateEnd = rec.mateEnd.getOrElse {
                  throw new IllegalArgumentException(s"Could not retrieve mate end for read '${rec.name}': Is the mate cigar present?")
                }
                def recOk  = rec.start <= strInterval.getStart && strInterval.getEnd <= rec.end
                def mateOk = rec.refIndex == rec.mateRefIndex && rec.mateStart <= strInterval.getStart && strInterval.getEnd <= mateEnd
                if (recOk || mateOk) {
                  true
                }
                else {
                  if (overlaps.nonEmpty) numNonSpaningReads += 1
                  false
                }
              }
            }
            else {
              detector.overlapsAny(new Interval(rec.refName, start, end))
            }
          }
        }
    }
  }
}

/** Metrics produced by `ReadGroupPerDuplexMolecularId` describing the counts per-strand for each UMI.
  *
  * @param mid the unique molecular identifier.
  * @param counts_a the counts for reads labelled `/A`.
  * @param counts_b the ocunts for reads labelled `/B`.
  * @param pass_filter true if the umi passed the min reads filter, false otherwise.
  */
case class UmiCounts(mid: String, counts_a: Long, counts_b: Long, pass_filter: Boolean) extends Metric
