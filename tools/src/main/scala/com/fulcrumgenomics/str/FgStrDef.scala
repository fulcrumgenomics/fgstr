package com.fulcrumgenomics.str

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.CommonsDef
import com.fulcrumgenomics.umi.ConsensusTags

object FgStrDef extends CommonsDef {

  /** Extract the molecular id from the record */
  def toMolecularId(rec: SamRecord, assignTag: String = ConsensusTags.MolecularId): Long = {
    val mi = rec[String](assignTag)
    val index = mi.lastIndexOf('/')
    require(index > 0, s"Read ${rec.name}'s $assignTag tag doesn't look like a duplex id: $mi")
    mi.substring(0, index).toLong
  }

  object DuplexFilters {
    /** Creates a set of duplex filters from a sequence of integers. Rules are as follows for `minReads`:
      * - The first value is always used across all read values
      * - The second value is used for AB filtering if provided, otherwise the first value
      * - The third value is used for BA filtering if provided, otherwise the second value, otherwise the first
      */
    def apply(minReads: Seq[Int]): DuplexFilters = {
      val filters = DuplexFilters(
        minReadsPerMolecule = minReads.head,
        minReadsAb          = minReads.take(2).last,
        minReadsBa          = minReads.last
      )
      require(filters.minReadsAb <= filters.minReadsPerMolecule, "min-reads values must be specified high to low.")
      require(filters.minReadsBa <= filters.minReadsAb, "min-reads values must be specified high to low.")
      filters
    }
  }

  /** Stores the minimum # of reads per-molecule, per-A-strand, and per-B-strand, where the minimum for the A-strand is
    * greater than the minimum for the B-strand.
    */
  case class DuplexFilters(minReadsPerMolecule: Int, minReadsAb: Int, minReadsBa: Int) {
    /** Returns the records to keep.  By default, filters out fragment reads. */
    def filter(recs: TraversableOnce[SamRecord], assignTag: String = ConsensusTags.MolecularId): Seq[SamRecord] = {
      val (pairs, _) = recs.toSeq.partition(_.paired)
      val r1s   = pairs.filter(_.firstOfPair)
      val total = r1s.length
      val (aStrand, bStrand) = {
        val a = r1s.count(r => r[String](assignTag).endsWith("/A"))
        val b = r1s.count(r => r[String](assignTag).endsWith("/B"))
        if (a > b) (a, b) else (b, a)
      }
      if (this.minReadsPerMolecule <= total && this.minReadsAb <= aStrand && this.minReadsBa <= bStrand) pairs.toSeq else Seq.empty
    }
  }
}
