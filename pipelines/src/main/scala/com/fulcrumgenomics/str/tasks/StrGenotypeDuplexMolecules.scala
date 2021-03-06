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

import com.fulcrumgenomics.FgBioDef.PathToBam
import dagr.tasks.DagrDef.{PathToFasta, PathToIntervals}

import scala.collection.mutable.ListBuffer

class StrGenotypeDuplexMolecules(input: PathToBam,
                                 output: PathToBam,
                                 ref: PathToFasta,
                                 intervals: PathToIntervals,
                                 perStrand: Boolean = false,
                                 minCumulativeFrequency: Option[Double] = None)
  extends FgStrTask {
  requires(1, "2g")

  protected def addFgStrArgs(buffer: ListBuffer[Any]): Unit = {
    buffer.append("-i", input)
    buffer.append("-o", output)
    buffer.append("-r", ref)
    buffer.append("-l", intervals)
    buffer.append("-s", perStrand)
    minCumulativeFrequency.foreach(m => buffer.append("-m", m))
  }
}