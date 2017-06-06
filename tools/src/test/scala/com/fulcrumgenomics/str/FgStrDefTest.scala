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

package com.fulcrumgenomics.str

import com.fulcrumgenomics.str.FgStrDef.DuplexFilters
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

class FgStrDefTest extends UnitSpec {

  "DuplexFilters" should "not reject read pairs if one end of the source molecule is not covered with with minReads=0" in {
    val builder = new SamBuilder
    builder.addPair(contig=0, start1=1, unmapped2=true, attrs=Map("MI" -> "1234/A"))
    DuplexFilters(minReads=Seq(0)).filter(builder.iterator.toSeq).length shouldBe 2
  }

  it should "reject read pairs if one end of the source molecule is not covered with minReads=1" in {
    val builder = new SamBuilder
    builder.addPair(contig=0, start1=1, unmapped2=true, attrs=Map("MI" -> "1234/A"))
    DuplexFilters(minReads=Seq(1)).filter(builder.iterator.toSeq) shouldBe 'empty
  }

  it should "support seeing reads only on one strand" in {
    val builder = new SamBuilder
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    // "--min-reads 1 1 0" allows there to be zero reads on the BA strand
    DuplexFilters(minReads=Seq(1, 1, 0)).filter(builder.iterator.toSeq).length shouldBe 2
  }

  it should "filter out fragment reads" in {
    val builder = new SamBuilder
    builder.addFrag(contig=0, start=1, attrs=Map("MI" -> "1234/A"))
    DuplexFilters(minReads=Seq(0)).filter(builder.iterator.toSeq) shouldBe 'empty
  }

  it should "support various incantations minReads" in {
    val builder = new SamBuilder

    // three total reads, two AB reads, one BA read
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, strand1=SamBuilder.Minus, strand2=SamBuilder.Plus, attrs=Map("MI" -> "1234/B"))

    // rejects since BA has only one read
    DuplexFilters(minReads=Seq(2)).filter(builder.iterator.toSeq) shouldBe 'empty

    // accepts since the last value (for AB) is repeated for BA
    DuplexFilters(minReads=Seq(2, 1)).filter(builder.iterator.toSeq).length shouldBe 6

    // accepts
    DuplexFilters(minReads=Seq(3, 2, 1)).filter(builder.iterator.toSeq).length shouldBe 6

    // rejects, not enough overall reads
    DuplexFilters(minReads=Seq(4, 2, 1)).filter(builder.iterator.toSeq) shouldBe 'empty
  }
}
