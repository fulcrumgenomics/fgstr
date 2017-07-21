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

import com.fulcrumgenomics.commons.CommonsDef.{DirPath, yieldAndThen}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.sopt.{arg, clp}
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.Pipeline
import dagr.tasks.DagrDef.PathToFasta

@clp(group=classOf[Pipelines], description=
  """
    |Multi-sample pipeline to genotype STRs using Duplex Sequencing.
  """)
class MultiSampleGenotypeWithDuplexSequencing
(
  @arg(flag='i', doc="Input directory with one sample per-directory") val input: DirPath,
  @arg(flag='o', doc="Output directory.") val output: DirPath,
  @arg(flag='r', doc="Path to the reference FASTA.") val ref: PathToFasta,
  @arg(flag='l', doc="Directory with interval lists, named <sample-name>.interval_list.") val intervals: DirPath,
  @arg(flag='A', doc="The read structure for read one.") val readStructureReadOne: String,
  @arg(flag='B', doc="The read structure for read two.") val readStructureReadTwo: String,
  @arg(flag='t', doc="Path to a temporary directory.  Use output if none is given.") val tmp: Option[DirPath] = None,
  @arg(flag='1', doc="Input fastq file suffix (optionally gzipped) for read 1.") val fastq1Suffix: String = "_L001_R1_001.fastq.gz",
  @arg(flag='2', doc="Input fastq file suffix (optionally gzipped) for read 2.") val fastq2Suffix: String = "_L001_R2_001.fastq.gz",
  @arg(flag='u', doc="The tag containing the raw UMI.")           val umiTag: String = "RX",
  @arg(flag='m', doc="Minimum mapping quality to include reads.") val minMapQ: Int = 10,
  @arg(flag='x', doc="The allowable number of edits between UMIs.") val edits: Int = 1,
  @arg(flag='M', minElements=1, maxElements=3, doc="The minimum number of raw reads per source molecule.")
  val minReads: Seq[Int] = Seq(1),
  @arg(flag='s', doc="Call genotypes per-duplex-strand") val perStrand: Boolean = false,
  @arg(          doc="Keep intermediate files when genotyping.") val keepIntermediates: Boolean = false

)  extends Pipeline(outputDirectory=Some(output.getParent)) {
  import MultiSampleGenotypeWithDuplexSequencing.SingletonOrException

  Io.assertListable(input)
  Io.assertReadable(ref)
  Io.assertListable(intervals)

  def build(): Unit = {
    Io.mkdirs(output)

    input.toFile.listFiles
      .filter(_.isDirectory)
      .foreach { inputSampleDir =>
        val sampleName      = inputSampleDir.getName
        val fastq1          = inputSampleDir.listFiles.filter(_.isFile).filter(_.getName.endsWith(fastq1Suffix)).toSeq.getSingleton
        val fastq2          = inputSampleDir.listFiles.filter(_.isFile).filter(_.getName.endsWith(fastq2Suffix)).toSeq.getSingleton
        val outputSampleDir = output.resolve(sampleName)
        val tmpSampleDir    = tmp.map(_.resolve(sampleName)).map { t => yieldAndThen(t)(Io.mkdirs(t)) }
        val sampleIntervals = intervals.resolve(sampleName + ".interval_list")

        Io.mkdirs(outputSampleDir)

        root ==> new GenotypeWithDuplexSequencing(
          fastq1               = List(fastq1.toPath),
          fastq2               = List(fastq2.toPath),
          output               = outputSampleDir.resolve(sampleName),
          tmp                  = tmpSampleDir,
          sample               = sampleName,
          library              = sampleName,
          platformUnit         = List(sampleName),
          readStructureReadOne = readStructureReadOne,
          readStructureReadTwo = readStructureReadTwo,
          ref                  = ref,
          intervals            = sampleIntervals,
          umiTag               = umiTag,
          minMapQ              = minMapQ,
          edits                = edits,
          minReads             = minReads,
          perStrand            = perStrand,
          keepIntermediates    = keepIntermediates
        )
      }
  }
}

object MultiSampleGenotypeWithDuplexSequencing {
  private implicit class SingletonOrException[T](things: Seq[T]) {
    def getSingleton: T = things match {
      case Seq() => throw new IllegalArgumentException("Did not find any elements")
      case Seq(thing) => thing
      case _ => throw new IllegalArgumentException(s"Found ${things.length} things, expected a single element")
    }
  }
}