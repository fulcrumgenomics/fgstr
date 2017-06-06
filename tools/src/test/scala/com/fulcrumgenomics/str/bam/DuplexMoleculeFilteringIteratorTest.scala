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

package com.fulcrumgenomics.str.bam

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

class DuplexMoleculeFilteringIteratorTest extends UnitSpec {

  def toIterator(iterator: Iterator[SamRecord], minReads: Seq[Int] = Seq(1)): DuplexMoleculeFilteringIterator =
    DuplexMoleculeFilteringIterator(iterator=iterator, minReads=minReads)

  "DuplexMoleculeFilteringIterator" should "throw an exception if the MI tag is null" in {
    val builder = new SamBuilder
    builder.addPair()
    an[Exception] should be thrownBy toIterator(builder.iterator).hasNext()
  }

  it should "throw an exception if the MI tag doesn't have a /suffix" in {
    val builder = new SamBuilder
    builder.addPair(attrs=Map("MI" -> "1234"))
    an[Exception] should be thrownBy toIterator(builder.iterator).hasNext()
  }

  it should "apply filters across a range of duplex filters" in {
    val builder = new SamBuilder

    // three total reads, two AB reads, one BA read
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, strand1=SamBuilder.Minus, strand2=SamBuilder.Plus, attrs=Map("MI" -> "1234/B"))

    // rejects since BA has only one read
    toIterator(builder.iterator, minReads=Seq(2)).isEmpty shouldBe true

    // accepts since the last value (for AB) is repeated for BA
    toIterator(builder.iterator, minReads=Seq(2, 1)).length shouldBe 6

    // accepts
    toIterator(builder.iterator, minReads=Seq(3, 2, 1)).length shouldBe 6

    // rejects, not enough overall reads
    toIterator(builder.iterator, minReads=Seq(4, 2, 1)).isEmpty shouldBe true
  }

  it should "support multiple duplex source molecules" in {
    val builder = new SamBuilder

    // fragment read
    builder.addFrag(unmapped=true, attrs=Map("MI" -> "1234/A"))

    // three total reads, two AB reads, one BA read (kept)
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, strand1=SamBuilder.Minus, strand2=SamBuilder.Plus, attrs=Map("MI" -> "1234/B"))

    // four total reads, three AB reads, one BA read (kept)
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "4321/A"))
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "4321/A"))
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "4321/A"))
    builder.addPair(contig=0, start1=1, start2=1, strand1=SamBuilder.Minus, strand2=SamBuilder.Plus, attrs=Map("MI" -> "4321/B"))

    // two total reads, one AB read and one BA read (filtered)
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "9876/A"))
    builder.addPair(contig=0, start1=1, start2=1, strand1=SamBuilder.Minus, strand2=SamBuilder.Plus, attrs=Map("MI" -> "9876/B"))

    val iter = toIterator(builder.iterator, minReads=Seq(3, 2, 1))
    val recs = iter.toIndexedSeq // force all reads to be consumed

    iter.isEmpty shouldBe true

    // # seen
    builder.iterator.length shouldBe 19
    iter.seen shouldBe 19
    // # kept
    recs.length shouldBe 14
    iter.kept shouldBe 14
    // # filtered
    iter.filtered shouldBe 5
  }
}
