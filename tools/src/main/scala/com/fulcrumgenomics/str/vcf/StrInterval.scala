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
package com.fulcrumgenomics.str.vcf

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.CommonsDef.PathToIntervals
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.str.vcf.StrInterval.StrAllele
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.variant.variantcontext.{Allele, VariantContext}

object StrInterval {
  object StrAllele extends LazyLogging {
    val NoCall: StrAllele = StrAllele(Allele.NO_CALL, 0, 0, 0)

    def toCalls(str: StrInterval, ctx: VariantContext, counts: Seq[Int], warn: Boolean = true): Seq[StrAllele] = {
      val alleles = ctx.getAlleles.toSeq
      require(alleles.length == counts.length, s"# of alleles '${alleles.length}' and # of counts '${counts.length}' did not match")
      val refAlleleLength = ctx.getReference.length()

      if (warn && str.unitLength * str.refLength != refAlleleLength) {
        logger.warning(s"Mismatch between reference repeat length in the interval list '${str.unitLength * str.refLength}' and vcf '$refAlleleLength' length")
      }

      alleles.zip(counts)
        .filter(_._2 > 0) // ignore zero counts
        .map { case (allele, count) =>
        // NB: the VCF's reference allele may not have the same length as the STR "reference" length, so we calculate the #
        // of bases different between the VCF's reference allele and this allele, then add that to the STR "reference"
        // length found in the interval list (unit-length times ref-length).
        val baseDiff     = allele.length() - refAlleleLength
        val alleleLength = baseDiff + (str.refLength * str.unitLength)
        StrAllele(allele=allele, alleleLength=alleleLength, count=count, str=str)
      }.sortBy(-_.count)
    }

    def apply(allele: Allele, alleleLength: Int, count: Int, str: StrInterval): StrAllele = {
      StrAllele(allele, alleleLength, count, str.unitLength)
    }
  }

  /** A single STR allele, including the number of repeats and count */
 case class StrAllele(allele: Allele, alleleLength: Int, count: Int, unitLength: Int) {
    def repeatLength: Double = {
      if (Allele.NO_CALL == allele) 0 else this.alleleLength / this.unitLength.toDouble
    }

    def toGenotype: String = {
      if (Allele.NO_CALL == allele) {
        "0"
      }
      else if (alleleLength % unitLength == 0) {
        f"${alleleLength / unitLength}%d"
      }
      else {
        f"$repeatLength%.2f"
      }
    }
  }

  /** Reads in the interval list with extra STR info. */
  def loadIntervals(path: PathToIntervals): Iterator[StrInterval] = {
    val intervals = IntervalList.fromFile(path.toFile)
    intervals.map { interval =>
      val strInfo = StrInterval(
        chrom      = interval.getContig,
        start      = interval.getStart,
        end        = interval.getEnd,
        // these will be set below
        unitLength = 0,
        refLength  = 0,
        name       = ""
      )

      interval.getName.split(",").toList match {
        case unitLength :: refLength :: name :: Nil =>
          strInfo.copy(unitLength=unitLength.toInt, refLength=refLength.toInt, name=name)
        case unitLength :: refLength :: name :: truthCalls =>
          strInfo.copy(unitLength=unitLength.toInt, refLength=refLength.toInt, name=name, truthCalls=truthCalls.map(_.toFloat))
        case _ =>
          throw new IllegalArgumentException(s"Interval name improperly formatted for interval: $interval")
      }
    }
  }
}

case class StrInterval(chrom: String,
                       start: Int,
                       end: Int,
                       unitLength: Int,
                       refLength: Int,
                       name: String,
                       truthCalls: Seq[Float] = Seq.empty)
  extends Interval(chrom, start, end, true, name) {
  override def toString: String  = productIterator.flatMap {
    case x: Option[_] => x
    case x: Seq[_]    => if (x.isEmpty) None else Some(x.mkString(","))
    case x            => Some(x)
  }.mkString("\t")

  /** Outputs a formatted string with the given set of (called) alleles */
  def toLongString(alleles: Seq[StrAllele]): String = {
    val genotypes = alleles.map { allele =>
      s"${allele.toGenotype}:${allele.count}"
    }.mkString(",")
    s"$this\t$genotypes"
  }
}