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

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SAMReadGroupRecord
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import htsjdk.samtools.util.{Interval, IntervalList, Iso8601Date}

class ReadGroupPerDuplexMolecularIdTest extends UnitSpec {

  private case class ReadGroupsAndRecords(readGroups: Seq[SAMReadGroupRecord], records: Seq[SamRecord])

  private def run(builder: SamBuilder, intervals: Option[PathToIntervals] = None, minReads:Seq[Int]= Seq(0), perStrand: Boolean=false): ReadGroupsAndRecords = {
    val input      = builder.toTempFile()
    val output     = makeTempFile("output", ".bam")
    executeFgstrTool(new ReadGroupPerDuplexMolecularId(input=input, output=output, intervals=intervals, minReads=minReads, perStrand=perStrand))
    val reader     = SamSource(output)
    val readGroups = reader.header.getReadGroups.toList
    val records    = reader.toIndexedSeq
    reader.close()
    ReadGroupsAndRecords(readGroups, records)
  }

  "ReadGroupPerDuplexMolecularId" should "fail if the input has no read groups" in {
    val builder = new SamBuilder()
    // remove the read group
    builder.header.setReadGroups(Iterator().toJavaList)
    val exception = intercept[IllegalArgumentException] { run(builder) }
    exception.getMessage should include ("Found no read groups")
  }

  it should "fail if the input has read groups with different sample/library combinations" in {
    val builder = new SamBuilder()
    def toRg(id: String): SAMReadGroupRecord = {
      val r = new SAMReadGroupRecord(id)
      r.setSample(s"Sample-for-$id")
      r
    }
    // Set the header to have two read groups
    val readGroups = Iterator(toRg("A"), toRg("B"))
    builder.header.setReadGroups(readGroups.toJavaList)
    val exception = intercept[IllegalArgumentException] { run(builder) }
    exception.getMessage should include ("Expected all read groups to have the same sample/library")
  }

  it should "fail if the molecular id is not an integer" in {
    val builder = new SamBuilder(readGroupId=Some("ID"))
    builder.addFrag(attrs=Map("MI" -> "ABCD/A"), unmapped=true)
    an[NumberFormatException] should be thrownBy run(builder)
  }

  it should "copy the input read group metadata over when one read group is present" in {
    val builder = new SamBuilder(readGroupId=Some("1234"))

    builder.rg.setSample("Sample")
    builder.rg.setLibrary("Library")
    builder.rg.setPlatformUnit("PlatformUnit")
    builder.rg.setPlatform("Platform")
    builder.rg.setRunDate(new Iso8601Date("2017-01-01T17:00:00-0700"))
    builder.rg.setSequencingCenter("SequencingCenter")
    builder.rg.setDescription("Description")
    builder.rg.setPredictedMedianInsertSize(42)
    builder.rg.setPlatformModel("PlatformModel")
    builder.header.setReadGroups(Iterator(builder.rg).toJavaList)

    builder.addFrag(attrs=Map("MI" -> "1234/A"), unmapped=true)

    val ReadGroupsAndRecords(readGroups, _) = run(builder)

    // update to make the tests below easier
    builder.rg.setSample("Sample-1234")
    builder.rg.setLibrary("Library-1234")

    readGroups.length shouldBe 1
    readGroups.head.equivalent(builder.rg) shouldBe true
  }

  private def toReadGroup(id: String, extra: String=""): SAMReadGroupRecord = {
    val rg = new SAMReadGroupRecord(id)

    rg.setSample("Sample")
    rg.setLibrary("Library")
    rg.setPlatformUnit("PlatformUnit" + extra)
    rg.setPlatform("Platform" + extra)
    rg.setRunDate(new Iso8601Date("2017-01-01T17:00:00-0700"))
    rg.setSequencingCenter("SequencingCenter" + extra)
    rg.setDescription("Description" + extra)
    rg.setPredictedMedianInsertSize(42)
    rg.setPlatformModel("PlatformModel" + extra)

    rg
  }

  it should "create an output read group when two read groups are present and have the same sample and library" in {
    val rgA = toReadGroup("A", extra="1")
    val rgB = toReadGroup("B", extra="2")

    val builder = new SamBuilder(readGroupId=Some("A"))

    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A", "RG" -> "A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A", "RG" -> "B"))

    // Set the read groups!
    builder.header.setReadGroups(Iterator(rgA, rgB).toJavaList)

    val ReadGroupsAndRecords(readGroups, records) = run(builder, minReads=Seq(0))

    readGroups.length shouldBe 1
    readGroups.map(_.getId) should contain theSameElementsInOrderAs Seq("mid-1")
    val rg = readGroups.head

    rg.getSample shouldBe "Sample-1"
    rg.getLibrary shouldBe "Library-1"
    rg.getPlatformUnit shouldBe "PlatformUnit1,PlatformUnit2"
    rg.getPlatform shouldBe "Platform1,Platform2"
    rg.getRunDate.toString.nonEmpty shouldBe true
    rg.getSequencingCenter shouldBe "SequencingCenter1,SequencingCenter2"
    rg.getDescription shouldBe "Description1,Description2"
    rg.getPredictedMedianInsertSize shouldBe 42
    rg.getPlatformModel shouldBe "PlatformModel1,PlatformModel2"

    records.length shouldBe 4
    records.map(_.attributes("RG")) should contain theSameElementsInOrderAs Seq("mid-1", "mid-1", "mid-1", "mid-1")
  }

  it should "fail if there were multiple sample/libraries" in {
    val rgA = toReadGroup("A", extra="1")
    val rgB = toReadGroup("B", extra="2")

    rgA.setSample("Sample-A")
    rgB.setSample("Sample-B")

    val builder = new SamBuilder(readGroupId=Some("A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A", "RG" -> "A"))

    // Set the read groups!
    builder.header.setReadGroups(Iterator(rgA, rgB).toJavaList)

    val exception = intercept[IllegalArgumentException] { run(builder, minReads=Seq(0)) }
    exception.getMessage should include ("Expected all read groups to have the same sample/library")
  }

  it should "output read groups per-strand" in {
    val builder = new SamBuilder(readGroupId=Some("ID"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/B"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "2/A"))

    val ReadGroupsAndRecords(readGroups, records) = run(builder, minReads=Seq(0), perStrand=true)

    readGroups.length shouldBe 3
    readGroups.map(_.getId) should contain theSameElementsAs Seq("mid-1/A", "mid-1/B", "mid-2/A")

    records.length shouldBe 8
    records.map(_.attributes("RG")) should contain theSameElementsInOrderAs Seq("mid-1/A", "mid-1/A", "mid-1/A", "mid-1/A", "mid-1/B", "mid-1/B", "mid-2/A", "mid-2/A")
  }

  it should "not care about the order of inputs" in {
    val builder = new SamBuilder(readGroupId=Some("ID"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "2/A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/B"))

    val ReadGroupsAndRecords(readGroups, records) = run(builder, minReads=Seq(0))

    readGroups.length shouldBe 2
    readGroups.head.getId shouldBe "mid-1"
    readGroups.last.getId shouldBe "mid-2"

    records.length shouldBe 6
    records.map(_.attributes("RG")) should contain theSameElementsInOrderAs Seq("mid-1", "mid-1", "mid-2", "mid-2", "mid-1", "mid-1")
  }

  it should "handle duplex reads" in {
    val builder = new SamBuilder(readGroupId=Some("ID"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/A"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "1/B"))
    builder.addPair(start1=1, start2=300, attrs=Map("MI" -> "2/A"))

    val ReadGroupsAndRecords(readGroups, records) = run(builder, minReads=Seq(0))

    readGroups.length shouldBe 2
    readGroups.head.getId shouldBe "mid-1"
    readGroups.last.getId shouldBe "mid-2"

    records.length shouldBe 6
    records.map(_.attributes("RG")) should contain theSameElementsInOrderAs Seq("mid-1", "mid-1", "mid-1", "mid-1", "mid-2", "mid-2")
  }

  it should "only output reads that overlap the given intervals" in {
    val builder = new SamBuilder(readGroupId=Some("ID"))

    val intervals = new IntervalList(builder.header)
    intervals.add(new Interval(builder.header.getSequence(0).getSequenceName, 1, 100))
    val ilist = makeTempFile("test.", ".interval_list")
    intervals.write(ilist.toFile)

    // fragment
    builder.addFrag(attrs=Map("MI" -> "1/A"), unmapped=true)
    builder.addFrag(attrs=Map("MI" -> "1/B"), unmapped=true)
    builder.addFrag(attrs=Map("MI" -> "2/A"), unmapped=true)

    // unmapped
    builder.addPair(attrs=Map("MI" -> "2/A"), unmapped1=true, unmapped2=true)
    builder.addPair(attrs=Map("MI" -> "2/B"), unmapped1=true, unmapped2=true)
    builder.addPair(attrs=Map("MI" -> "3/A"), unmapped1=true, unmapped2=true)

    // mapped in the intervals
    builder.addPair(attrs=Map("MI" -> "3/A"), contig=0, start1=1, start2=1)
    builder.addPair(attrs=Map("MI" -> "3/B"), contig=0, start1=100, start2=100)

    // mapped outside
    builder.addPair(attrs=Map("MI" -> "4/A"), contig=0, start1=200, start2=200)
    builder.addPair(attrs=Map("MI" -> "5/A"), contig=1, start1=100, start2=100)

    val ReadGroupsAndRecords(readGroups, records) = run(builder, minReads=Seq(0), intervals=Some(ilist))

    readGroups.length shouldBe 1
    readGroups.head.getId shouldBe "mid-3"

    records.length shouldBe 4
    records.foreach { record =>
      record.attributes("RG") shouldBe "mid-3"
    }
  }

  it should "suport --min-reads" in {
    val builder = new SamBuilder

    // three total reads, two AB reads, one BA read
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, attrs=Map("MI" -> "1234/A"))
    builder.addPair(contig=0, start1=1, start2=1, strand1=SamBuilder.Minus, strand2=SamBuilder.Plus, attrs=Map("MI" -> "1234/B"))

    // rejects since BA has only one read
    run(builder,  minReads=Seq(2)).records.length shouldBe 0

    // accepts since the last value (for AB) is repeated for BA
    run(builder,  minReads=Seq(2, 1)).records.length shouldBe 6

    // accepts
    run(builder,  minReads=Seq(3, 2, 1)).records.length shouldBe 6

    // rejects, not enough overall reads
    run(builder,  minReads=Seq(4, 2, 1)).records.length shouldBe 0
  }

  it should "still output a read group if the input BAM is empty" in {
    Seq(true, false).foreach { perStrand =>
      val builder = new SamBuilder
      val ReadGroupsAndRecords(readGroups, records) = run(builder, minReads=Seq(0), perStrand=perStrand)

      readGroups.length shouldBe 1
      records shouldBe 'empty

      // make sure they are the same!
      if (perStrand) {
        builder.rg.setSample("Sample-0/A")
      }
      else {
        builder.rg.setSample("Sample-0")
      }
      readGroups.head.equivalent(builder.rg) shouldBe true
    }
  }
}
