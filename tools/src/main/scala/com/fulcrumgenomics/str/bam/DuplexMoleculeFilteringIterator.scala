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
import com.fulcrumgenomics.str.FgStrDef.{DuplexFilters, toMolecularId}
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.FgBioDef._
import scala.collection.mutable

object DuplexMoleculeFilteringIterator {
  def apply(iterator: Iterator[SamRecord],
            minReads: Seq[Int],
            assignTag: String = ConsensusTags.MolecularId): DuplexMoleculeFilteringIterator = {
    new DuplexMoleculeFilteringIterator(iterator, DuplexFilters(minReads), assignTag)
  }
}

/** Filters records from the same source duplex molecule based on coverage of the molecule, A-strand, and B-strand
  * respectively. */
class DuplexMoleculeFilteringIterator(iterator: Iterator[SamRecord],
                                      filters: DuplexFilters,
                                      assignTag: String = ConsensusTags.MolecularId
                                     ) extends Iterator[SamRecord] {
  private val bufferedIterator = iterator.bufferBetter
  private val nextRecords      = mutable.Queue[SamRecord]()
  private var _seen            = 0L
  private var _kept            = 0L

  // Get the first batch of reads
  this.advance()

  /** Returns true if there are more records, false otherwise. */
  override def hasNext(): Boolean = {
    if (this.nextRecords.isEmpty && this.bufferedIterator.hasNext) this.advance()
    this.nextRecords.nonEmpty
  }

  /** Returns the next record not filtered. */
  override def next(): SamRecord = {
    if (!hasNext) throw new NoSuchElementException("Calling next() when hasNext() is false")
    nextRecords.dequeue()
  }

  /** The number of records seen. */
  def seen: Long = this._seen

  /** The number of records kept. */
  def kept: Long = this._kept

  /** The number of records filtered. */
  def filtered: Long = seen - kept

  /** Tries to get the next set of reads from a source molecule that isn't filtere.d */
  private def advance(): Unit = {
    require(this.nextRecords.isEmpty, "Calling advance() when nextRecords is not empty")
    while (this.bufferedIterator.hasNext && this.nextRecords.isEmpty) {
      val recs           = nextGroupOfRecords()
      val recsToKeep     = filters.filter(recs)
      this.nextRecords ++= recsToKeep
      this._kept        += recsToKeep.length
      this._seen        += recs.length
    }
  }

  /**
    * Consumes the next group of records from the input iterator, based on molecule id
    * and returns them as a Seq.
    */
  private def nextGroupOfRecords(): Seq[SamRecord] = {
    val idToMatch = toMolecularId(this.bufferedIterator.head, this.assignTag)
    this.bufferedIterator.takeWhile(toMolecularId(_) == idToMatch).toSeq
  }
}