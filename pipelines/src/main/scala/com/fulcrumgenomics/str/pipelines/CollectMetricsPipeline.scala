package com.fulcrumgenomics.str.pipelines

import dagr.commons.io.Io
import dagr.core.cmdline.Pipelines
import dagr.core.tasksystem.Pipeline
import dagr.sopt.{arg, clp}
import dagr.tasks.DagrDef.{PathPrefix, PathToBam, PathToFasta, PathToIntervals}
import dagr.tasks.picard.{CollectAlignmentSummaryMetrics, CollectHsMetrics, CollectSequencingArtifactMetrics}

@clp(
  description =
    """
      |Pipeline to compute metrics.
      |
      |Will run CollectAlignmentSummaryMetrics, CollectHsMetrics, and CollectSequencingArtifactMetrics.
    """,
  group = classOf[Pipelines]
)
class CollectMetricsPipeline
(
  @arg(flag="o", doc="Input mapped BAM).") val input: PathToBam,
  @arg(flag="r", doc="Path to the reference FASTA.")              val ref: PathToFasta,
  @arg(flag="l", doc="Regions to analyze.")                       val intervals: PathToIntervals,
  @arg(flag="o", doc="Path prefix for output files.")             val output: PathPrefix
)
  extends Pipeline(Some(output.getParent)) {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output, parentMustExist=false)

  def build(): Unit = {
    root ==> new CollectAlignmentSummaryMetrics(in=input, ref=ref)
    root ==> new CollectHsMetrics(in=input, ref=ref, targets=intervals)
    root ==> new CollectSequencingArtifactMetrics(in=input, ref=ref, intervals=Some(intervals), minBq=Some(30), contextSize=Some(0))
  }
}

