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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.ClpGroups
import com.fulcrumgenomics.commons.CommonsDef.IteratorToJavaCollectionsAdapter
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.str.FgStrDef
import com.fulcrumgenomics.str.FgStrDef.DuplexFilters
import com.fulcrumgenomics.str.cmdline.FgStrTool
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector}
import htsjdk.samtools.{SAMFileHeader, SAMReadGroupRecord}

import scala.collection.mutable

@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Creates a BAM file with a read group per duplex molecule ID.
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
  @arg(flag='s', doc="Create a read group per-strand of a duplex molecule") val perStrand: Boolean = false
) extends FgStrTool with LazyLogging {

  Io.assertReadable(input)
  intervals.foreach(Io.assertReadable)
  Io.assertCanWriteFile(output)

  private val toMolecularId: SamRecord => String = rec => if (this.perStrand) rec[String](assignTag) else FgStrDef.toMolecularId(rec, assignTag).toString
  private val filters = DuplexFilters(minReads)

  validate(input.toAbsolutePath != Io.StdIn, "Reading from standard input is not supported")

  override def execute(): Unit = {
    // Pass 1: find all the molecular identifiers to use
    logger.info("Extracting molecular ids")
    val mids = {
      val progress = ProgressLogger(this.logger, unit=5e5.toInt, verb="Extracted molecule ids")
      val in       = SamSource(input, ref=ref)
      val iterator = toIterator(in)
      val ids      = new mutable.HashSet[String]()
      iterator.foreach { rec =>
        ids += toMolecularId(rec)
        progress.record()
        progress.record(rec)
        if (progress.getCount % 5e5.toInt == 0) {
          logger.info("# of molecules so far: " + ids.size)
        }
      }
      in.safelyClose()
      logger.info(s"Extracted ${ids.size} molecular ids across ${progress.getCount} raw reads")
      ids.toSeq.sorted
    }

    // Pass 2: update the header and records
    {
      logger.info("Writing to output")
      val in       = SamSource(input, ref=ref)
      val iterator = new DuplexMoleculeFilteringIterator(toIterator(in), filters, assignTag)
      val header   = toOutputHeader(in, mids)
      val out      = SamWriter(output, header, ref=ref)

      if (mids.nonEmpty) {
        iterator.foreach { rec =>
          rec("RG") = s"mid-${toMolecularId(rec)}"
          out += rec
        }
      }
      in.safelyClose()
      out.close()
    }
  }

  /** Create the output header. */
  private def toOutputHeader(in: SamSource, mids: Seq[String]): SAMFileHeader = {
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

    def toReadGroup(mid: String): SAMReadGroupRecord = {
      val readGroup = new SAMReadGroupRecord(s"mid-$mid")
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
      mids.map { mid => toReadGroup(mid) }
    }
    else {
      // We may not have ANY mids if there were no reads, or reads did not overlap the intervals,
      // so just create a "dummy" read group.
      Seq(toReadGroup(if (perStrand) "0/A" else "0"))
    }

    header.setReadGroups(readGroups.toIterator.toJavaList)
    header
  }

  /** Build the iterator we'll use based on whether or not we're restricting to a set of intervals */
  private def toIterator(in: SamSource): Iterator[SamRecord] = {
    intervals match {
      case None       => in.view.toIterator
      case Some(path) =>
        val ilist    = IntervalList.fromFile(path.toFile).uniqued(false)
        val detector = new OverlapDetector[Interval](0,0)
        detector.addAll(ilist.getIntervals, ilist.getIntervals)
        in.iterator.filter { rec: SamRecord =>
          if (!rec.paired || rec.unmapped || rec.mateUnmapped) false
          else {
            val (start, end) = if (rec.refIndex == rec.mateRefIndex) Bams.insertCoordinates(rec) else (rec.start, rec.end)
            detector.overlapsAny(new Interval(rec.refName, start, end))
          }
        }
    }
  }
}