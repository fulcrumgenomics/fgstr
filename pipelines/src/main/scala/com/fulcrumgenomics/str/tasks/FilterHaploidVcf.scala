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

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToVcf}
import dagr.api.models.{Cores, Memory}
import dagr.core.tasksystem.{FixedResources, PipeOut, ProcessTask}
import dagr.tasks.DataTypes.Vcf

import scala.collection.mutable.ListBuffer

class FilterHaploidVcf(input: PathToVcf,
                       stats: Option[FilePath] = None,
                       minQuality: Option[Double] = None,
                       spanningReadsRequired: Boolean = true,
                       maxLocusDepth: Option[Int] = Some(100000000)) // effectively infinite
  extends ProcessTask with FixedResources with PipeOut[Vcf] {

  requires(Cores(1), Memory("2g"))

  // Find it when building so we fail earlier
  private val script = HipStr.findHipStrScript(scriptName="filter_haploid_vcf.py")

  override def args: Seq[Any] = {
    val buffer = ListBuffer[Any]()

    buffer.append("python")
    buffer.append(script)
    buffer.append("--vcf", input)
    stats.foreach(s => buffer.append("--stats", s))
    minQuality.foreach(q => buffer.append("--min-call-qual", q))
    if (spanningReadsRequired) buffer.append("--no-spanning")
    maxLocusDepth.foreach(d => buffer.append("--max-loc-depth", d))

    buffer
  }
}