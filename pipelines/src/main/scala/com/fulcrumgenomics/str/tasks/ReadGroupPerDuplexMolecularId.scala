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

import com.fulcrumgenomics.FgBioDef.{PathToBam, PathToFasta}
import com.fulcrumgenomics.bam.api.SamOrder
import dagr.tasks.DagrDef.PathToIntervals

import scala.collection.mutable.ListBuffer

class ReadGroupPerDuplexMolecularId(input: PathToBam,
                                    output: PathToBam,
                                    ref: Option[PathToFasta] = None,
                                    intervals: Option[PathToIntervals] = None,
                                    minReads: Seq[Int] = Seq(1),
                                    perStrand: Boolean = false,
                                    samOrder: Option[SamOrder] = Some(SamOrder.Coordinate),
                                    span: Boolean = false)
  extends FgStrTask {

  protected def addFgStrArgs(buffer: ListBuffer[Any]): Unit = {
    buffer.append("-i", input)
    buffer.append("-o", output)
    ref.foreach(r => buffer.append("-r", r))
    intervals.foreach(l => buffer.append("-l", l))
    buffer.append("-M"); buffer.append(minReads:_*)
    buffer.append("-s", perStrand)
    samOrder.foreach(s => buffer.append("-S", s))
    buffer.append("--span", span)
  }
}
