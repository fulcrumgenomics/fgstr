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

package com.fulcrumgenomics.str.cmdline

import com.fulcrumgenomics.cmdline.FgBioMain
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.commons.CommonsDef.unreachable
import com.fulcrumgenomics.sopt.Sopt

/**
  * Main program for FgStrTools that loads everything up and runs the appropriate sub-command
  */
object FgStrMain {
  /** The main method */
  def main(args: Array[String]): Unit = new FgStrMain().makeItSoAndExit(args)
}

class FgStrCommonArgs

class FgStrMain extends FgBioMain {
  override protected def name: String = "fgstr"

  /** The packages we wish to include in our command line **/
  override protected def packageList: List[String] = List[String]("com.fulcrumgenomics.str")

  /** A main method that returns an exit code instead of exiting. */
  override def makeItSo(args: Array[String]): Int = {
    // Turn down HTSJDK logging
    htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)

    val startTime = System.currentTimeMillis()
    val exit      = Sopt.parseCommandAndSubCommand[FgStrCommonArgs,FgStrTool](name, args, Sopt.find[FgStrTool](packageList)) match {
      case Sopt.Failure(usage) =>
        System.err.print(usage())
        1
      case Sopt.CommandSuccess(cmd) =>
        unreachable("CommandSuccess should never be returned by parseCommandAndSubCommand.")
      case Sopt.SubcommandSuccess(_, subcommand) =>
        val name = subcommand.getClass.getSimpleName
        try {
          printStartupLines(name, args)
          subcommand.execute()
          printEndingLines(startTime, name, true)
          0
        }
        catch {
          case ex: FailureException =>
            val banner = "#" * ex.message.map(_.length).getOrElse(80)
            logger.fatal(banner)
            logger.fatal("Execution failed!")
            ex.message.foreach(logger.fatal)
            logger.fatal(banner)
            printEndingLines(startTime, name, false)
            ex.exit
          case ex: Throwable =>
            printEndingLines(startTime, name, false)
            throw ex
        }
    }

    exit
  }

}
