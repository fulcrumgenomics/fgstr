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

import java.nio.file.Path

import dagr.core.execsystem.{Cores, Memory}
import dagr.tasks.picard.PicardTask
import htsjdk.samtools.ValidationStringency

import scala.collection.mutable.ListBuffer

object FgStrTask {
  /** Override to return the path to the fgstr-tools jar. */
  val FgStrJarConfigPath = "fgstr-tools.jar"
}

/**
  * Abstract base task that should be extended for any tool in the fgstr-tools toolset.
  */
abstract class FgStrTask(useAdvancedGcOptions: Boolean = true,
                         validationStringency: Option[ValidationStringency] = Some(ValidationStringency.LENIENT),
                         useAsyncIo: Boolean = false,
                         compressionLevel: Option[Int] = None,
                         createIndex: Option[Boolean] = Some(true),
                         createMd5File: Option[Boolean] = None)
  extends PicardTask(
    useAdvancedGcOptions = useAdvancedGcOptions,
    validationStringency = validationStringency,
    useAsyncIo           = useAsyncIo,
    compressionLevel     = compressionLevel,
    createIndex          = createIndex,
    createMd5File        = createMd5File) {

  // Default resource requirements
  requires(Cores(1), Memory("4G"))

  /** Override to return the path to the fgstr-tools jar. */
  override def jarPath: Path = configure[Path](FgStrTask.FgStrJarConfigPath)

  /** Override to "hide" the addPicardArgs being addFgStrArgs. */
  final override protected def addPicardArgs(buffer: ListBuffer[Any]): Unit = addFgStrArgs(buffer)

  /** To be implemented by sub-classes to add arguments to the command line. */
  protected def addFgStrArgs(buffer: ListBuffer[Any]): Unit
}
