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

import java.nio.file.{Files, Path}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.str.tasks.FgStrTask
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.commons.io.Io
import dagr.core.cmdline.Pipelines
import dagr.core.config.Configuration
import dagr.core.tasksystem.{Pipeline, ProcessTask}
import dagr.tasks.DagrDef.{PathPrefix, PathToFastq}
import dagr.tasks.DataTypes.SamOrBam
import dagr.tasks.fgbio.{ExtractUmisFromBam, FgBioTask}
import dagr.tasks.misc.MoveFile
import dagr.tasks.picard._
import htsjdk.samtools.SAMFileHeader.SortOrder

import scala.collection.mutable.ListBuffer

object PrepareUnmappedBamPipeline extends Configuration {
  private def requireConfigure(path: String): Unit = require(configure[Path](path) != null, s"Configuration path not found: $path")
}

@clp(group=classOf[Pipelines], description=
  """
    |Pipeline to prepare an unmapped BAM from one or more pair of FASTQs.
    |
    |The resulting unmapped BAM can be used as input to the `MapAndGroupRawReads` pipeline.
  """)
class PrepareUnmappedBamPipeline
(
  @arg(flag='1', doc="Input fastq file (optionally gzipped) for read 1.") val fastq1: List[PathToFastq],
  @arg(flag='2', doc="Input fastq file (optionally gzipped) for read 2.") val fastq2: List[PathToFastq],
  @arg(flag='o', doc="Output file prefix (e.g.dir/sample_name).")         val output: PathPrefix,
  @arg(flag='S', doc="The name of the sample.")                           val sample: String,
  @arg(flag='L', doc="The name of the library.")                          val library: String,
  @arg(flag='P', doc="The platform unit (@RG.PU).  Either one value or one per pair of FASTQs.") val platformUnit: List[String],
  @arg(flag='t', doc="Path to a temporary directory.  Use output if none is given.") val tmp: Option[DirPath] = None,
  @arg(flag='A', doc="The read structure for read one.") val readStructureReadOne: String,
  @arg(flag='B', doc="The read structure for read two.") val readStructureReadTwo: String
)  extends Pipeline(outputDirectory=Some(output.getParent), suffix=Some("." + library)) {
  import PrepareUnmappedBamPipeline.requireConfigure

  requireConfigure(FgBioTask.FgBioJarConfigPath)
  requireConfigure(PicardTask.PicardJarConfigPath)
  requireConfigure(FgStrTask.FgStrJarConfigPath)
  require(fastq1.length == fastq2.length, "--fastq1 and --fastq2 must have the same # of arguments")
  require(platformUnit.length == fastq1.length || platformUnit.length == 1, "--platform-unit must have a single value or the same # of values as --fastq1")

  private val prefix: String = output.getFileName.toString
  private val out = output.getParent
  val unmappedBamFile: PathPrefix = out.resolve(prefix + ".unmapped.bam")

  override def build(): Unit = {
    val tmpDir = tmp.getOrElse(out.resolve("tmp"))

    Io.assertReadable(fastq1 ++ fastq2)
    Io.assertCanWriteFile(out.resolve(prefix))
    Io.mkdirs(tmpDir)

    val platformUnits       = if (1 < platformUnit.length) platformUnit else List.tabulate(fastq1.length)(x => platformUnit.head)
    val inputs              = (fastq1, fastq2, platformUnits).zipped
    val toUnmappedBams      = ListBuffer[ProcessTask]()
    val unmappedBams        = ListBuffer[PathToBam]()

    ///////////////////////////////////////////////////////////////////////
    // Generate an unmapped BAM for each pair of fastqs input
    ///////////////////////////////////////////////////////////////////////
    inputs.foreach((fq1, fq2, pu) => {
      val bam = Files.createTempFile(tmpDir, "unmapped.", ".bam")
      unmappedBams += bam

      val metricsPrefix = if (inputs.size > 1) out.resolve(prefix + s".${unmappedBams.length}") else out.resolve(prefix)
      val fastqToSam    = new FastqToSam(fastq1=fq1, fastq2=Some(fq2), out=Io.StdOut, sample=sample, library=Some(library), readGroupName=None, platformUnit=Some(pu))
      val buffer1       = new FifoBuffer[SamOrBam]
      val extractUmis   = new ExtractUmisFromBam(
        in                = Io.StdIn,
        out               = Io.StdOut,
        readStructures    = Seq(readStructureReadOne, readStructureReadTwo),
        umiTags           = Seq(ConsensusTags.UmiBases),
        annotateNames     = false
      )
      val buffer2       = new FifoBuffer[SamOrBam]
      val markAdapters  = new MarkIlluminaAdapters(in=Io.StdIn, out=bam, prefix=Some(metricsPrefix))
      val toUnmappedBam = (fastqToSam | buffer1 | extractUmis | buffer2 | markAdapters).withName("FastqToUnmappedBam")
      toUnmappedBams += toUnmappedBam
      root ==> toUnmappedBam
    })

    ///////////////////////////////////////////////////////////////////////
    // Then either merge all the input BAMs, or if we have just a single
    // one, then just rename it
    ///////////////////////////////////////////////////////////////////////
    unmappedBams.toList match {
      case uBam :: Nil =>
        toUnmappedBams.head ==> new MoveFile(uBam, unmappedBamFile).withName("Move unmapped BAM")
      case _ =>
        val mergeUnmappedSams = new MergeSamFiles(in=unmappedBams, out=unmappedBamFile, sortOrder=SortOrder.queryname)
        toUnmappedBams.foreach(_ ==> mergeUnmappedSams)
        unmappedBams.foreach(b => mergeUnmappedSams ==> new DeleteBam(b))
    }
  }
}

