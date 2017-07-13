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

import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.str.vcf.StrInterval.StrAllele
import com.fulcrumgenomics.testing.{UnitSpec, VariantContextSetBuilder}
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.variant.variantcontext.{Allele, Genotype, VariantContext}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.ReferenceSetBuilder

class StrGenotypeDuplexMoleculesTest extends UnitSpec {

  private val ref: PathToFasta = {
    val dict = new VariantContextSetBuilder().header.getSequenceDictionary
    val builder = new ReferenceSetBuilder
    dict.getSequences.foreach { sequence =>
      builder.add(sequence.getSequenceName).add("A", 1000) // NB: relies on the sequence dictionaries not being compared
    }
    builder.toTempFile()
  }

  "StrAllele.toGenotype" should "format as an integer if the # of repeat units is a multiple of the unit length" in {
    StrAllele(null, 10, 1, 10).toGenotype shouldBe "1"
    StrAllele(null, 10, 1, 5).toGenotype  shouldBe "2"
    StrAllele(null, 9,  1, 3).toGenotype  shouldBe "3"
  }

  it should "format as an float if the # of repeat units is not a multiple of the unit length" in {
    StrAllele(null, 10, 1, 3).toGenotype shouldBe "3.33"
    StrAllele(null, 10, 1, 4).toGenotype shouldBe "2.50"
  }

  private object Outputs {
    def apply(output: PathPrefix): Outputs = {
      val vcf = PathUtil.pathTo(output + ".vcf.gz")
      new Outputs(
        variants = readVcfRecs(vcf).toIndexedSeq,
        vcf      = vcf,
        info     = PathUtil.pathTo(output + ".info.txt"),
        pdf      = PathUtil.pathTo(output + ".pdf")
      )
    }
  }
  private case class Outputs(variants: IndexedSeq[VariantContext], vcf: PathToVcf, info: FilePath, pdf: FilePath)

  private val intervals: PathToIntervals = {
    val dict = new VariantContextSetBuilder().header.getSequenceDictionary
    val ilist = makeTempFile("test.", ".interval_list")
    val intvs = new IntervalList(dict)
    intvs.add(new Interval(dict.getSequence(0).getSequenceName, 1, 3, false, "3,1,NAME1"))
    intvs.write(ilist.toFile)
    ilist
  }

  private def getStrGenotype(genotype: Genotype): Seq[Double] = {
    genotype.hasAnyAttribute("STR_GT") shouldBe true
    genotype.getAttributeAsString("STR_GT", ".").split(',').map(_.toDouble)
  }

  private def run(builder: VariantContextSetBuilder,
                  intervals: PathToIntervals=this.intervals,
                  ref: PathToFasta=this.ref,
                  minCumulativeFrequency: Double = 0.9,
                  perStrand: Boolean = false,
                  skipPlots: Boolean = true): Outputs = {
    val input  = builder.toTempFile()
    val output = makeTempFile("test.", ".prefix")
    val tool   = new StrGenotypeDuplexMolecules(input=input, output=output, intervals=intervals, ref=ref,
      minCumulativeFrequency=minCumulativeFrequency, perStrand=perStrand, skipPlots=skipPlots)
    val logging = executeFgstrTool(tool)
    Outputs(output)
  }

  "StrGenotypeDuplexMolecules" should "fail if the interval list's name field is in the wrong format" in {
    val builder   = new VariantContextSetBuilder(sampleNames=Seq("mid-1234"))
    val dict      = builder.header.getSequenceDictionary
    val ilist     = makeTempFile("test.", ".interval_list")

    val intervals = new IntervalList(dict)
    intervals.add(new Interval(dict.getSequence(0).getSequenceName, 100, 200))
    intervals.add(new Interval(dict.getSequence(0).getSequenceName, 100, 200, false, "4,NAME"))
    intervals.write(ilist.toFile)

    val exception = intercept[IllegalArgumentException] { run(builder, ilist) }
    exception.getMessage should include ("Interval name improperly formatted")
  }

  it should "fail if the sample names do not have '-<mid>' suffixes" in {
    val names   = List("mid-1234", "mid", "mid-2345")
    val builder = new VariantContextSetBuilder(sampleNames=names)

    names.foreach { name =>
      builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAAAAA"), sampleName=Some(name))
    }

    an[Exception] should be thrownBy run(builder, intervals)
  }

  it should "skip variants that do not overlap the input intervals" in {
    val builder   = new VariantContextSetBuilder(sampleNames=Seq("mid-1234"))
    val dict      = builder.header.getSequenceDictionary
    val ilist     = makeTempFile("test.", ".interval_list")

    val intervals = new IntervalList(dict)
    intervals.add(new Interval(dict.getSequence(0).getSequenceName, 100, 200, false, "4,10,NAME"))
    intervals.write(ilist.toFile)

    // FIXME: remove genotypeAlleles=... when this PR is merged: https://github.com/fulcrumgenomics/fgbio/pull/245
    builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA"), genotypeAlleles=List(Allele.NO_CALL_STRING)) // before the interval
    builder.addVariant(refIdx=0, start=201, variantAlleles=List("AAA", "AAAA"), genotypeAlleles=List(Allele.NO_CALL_STRING)) // past the interval
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("AAA", "AAAA"), genotypeAlleles=List(Allele.NO_CALL_STRING)) // next chromosome

    val outputs  = run(builder, ilist)
    val variants = outputs.variants

    variants.length shouldBe 1
    variants.count(_.getNoCallCount == 0) shouldBe 0
  }

  it should "output two variants when overlapping STR intervals are given" in {
    val builder   = new VariantContextSetBuilder(sampleNames=List("mid-1234"))

    builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAAAAA"))

    val intervals = {
      val dict = builder.header.getSequenceDictionary
      val ilist = makeTempFile("test.", ".interval_list")
      val intvs = new IntervalList(dict)
      intvs.add(new Interval(dict.getSequence(0).getSequenceName, 1, 3, false, "3,1,NAME1"))
      intvs.add(new Interval(dict.getSequence(0).getSequenceName, 1, 6, false, "3,2,NAME2"))
      intvs.write(ilist.toFile)
      ilist
    }

    val outputs  = run(builder, intervals)
    val variants = outputs.variants

    variants.length shouldBe 2
    variants.count(_.getNoCallCount == 0) shouldBe 2
  }


  // NB: also tests creating plots
  it should "produce various heterozygous and homozygous calls" in {
    val builder   = new VariantContextSetBuilder(sampleNames=List("S-1", "S-2"))

    // homozygous reference
    builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA"), genotypeAlleles=List("AAA"), sampleName=Some("S-1"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA"), genotypeAlleles=List("AAA"), sampleName=Some("S-2"))

    // heterozygous reference
    builder.addVariant(refIdx=1, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAA"), sampleName=Some("S-1"))
    builder.addVariant(refIdx=1, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAAAAA"), sampleName=Some("S-2"))

    // homozygous alternate
    builder.addVariant(refIdx=2, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAAAAA"), sampleName=Some("S-1"))
    builder.addVariant(refIdx=2, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAAAAA"), sampleName=Some("S-2"))

    // heterozygous alternate
    builder.addVariant(refIdx=3, start=1, variantAlleles=List("AAA", "AAAAAA", "AAAAAAAAA"), genotypeAlleles=List("AAAAAA"), sampleName=Some("S-1"))
    builder.addVariant(refIdx=3, start=1, variantAlleles=List("AAA", "AAAAAA", "AAAAAAAAA"), genotypeAlleles=List("AAAAAAAAA"), sampleName=Some("S-2"))

    val dict = new VariantContextSetBuilder().header.getSequenceDictionary
    val ilist = makeTempFile("test.", ".interval_list")
    val intvs = new IntervalList(dict)
    intvs.add(new Interval(dict.getSequence(0).getSequenceName, 1, 3, false, "3,1,NAME1"))
    intvs.add(new Interval(dict.getSequence(1).getSequenceName, 1, 6, false, "3,1,NAME2"))
    intvs.add(new Interval(dict.getSequence(2).getSequenceName, 1, 6, false, "3,1,NAME3"))
    intvs.add(new Interval(dict.getSequence(3).getSequenceName, 1, 9, false, "3,1,NAME4"))
    intvs.write(ilist.toFile)

    val outputs  = run(builder, ilist, skipPlots=false)
    val variants = outputs.variants

    variants.length shouldBe 4

    val validateMethods: Seq[Genotype => Boolean] = Seq(
      g => g.isHomRef,
      g => g.isHet && !g.isHetNonRef,
      g => g.isHomVar,
      g => g.isHetNonRef
    )

    val validateStrGenotypes: Seq[Genotype => Unit] = Seq(
      g => getStrGenotype(g) should contain theSameElementsInOrderAs Seq(1.0, 1.0),
      g => getStrGenotype(g) should contain theSameElementsInOrderAs Seq(1.0, 2.0),
      g => getStrGenotype(g) should contain theSameElementsInOrderAs Seq(2.0, 2.0),
      g => getStrGenotype(g) should contain theSameElementsInOrderAs Seq(2.0, 3.0)
    )

    variants.zip(validateMethods).foreach { case (ctx, validateMethod) =>
      ctx.getGenotypes.length shouldBe 1
      validateMethod(ctx.getGenotype(0)) shouldBe true
    }

    variants.zip(validateStrGenotypes).foreach { case (ctx, validateStrGenotype) =>
      validateStrGenotype(ctx.getGenotype(0))
    }

    Io.assertReadable(outputs.info)
    Io.assertReadable(outputs.pdf)

    val actualLines = Io.readLines(outputs.info). toSeq
    val expectedLines = Seq(
      "chr1\t1\t3\t3\t1\tNAME1\t1:2",     // two occurrences of 1 unit
      "chr2\t1\t6\t3\t1\tNAME2\t1:1,2:1", // one occurrence of 1 unit, one occurrence of 2 units
      "chr3\t1\t6\t3\t1\tNAME3\t2:2",     // two occurrences of 2 unit
      "chr4\t1\t9\t3\t1\tNAME4\t2:1,3:1"  // one occurrence of 2 units, one occurrence of 3 units
    )
    actualLines should contain theSameElementsInOrderAs expectedLines
  }

   it should "handle an interval list with known calls" in {
    val builder   = new VariantContextSetBuilder(sampleNames=List("mid-1", "mid-2"))

    builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAA"), sampleName=Some("mid-1"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=List("AAA", "AAAAAA"), genotypeAlleles=List("AAAAAA"), sampleName=Some("mid-2"))

    val intervals = {
      val dict = builder.header.getSequenceDictionary
      val ilist = makeTempFile("test.", ".interval_list")
      val intvs = new IntervalList(dict)
      intvs.add(new Interval(dict.getSequence(0).getSequenceName, 1, 3, false, "3,1,NAME1,1,2"))
      intvs.write(ilist.toFile)
      ilist
    }

    val outputs  = run(builder, intervals)
    val variants = outputs.variants

    variants.length shouldBe 1
    variants.head.getGenotypes.length shouldBe 1
    val genotype = variants.head.getGenotype(0)
    genotype.isHet shouldBe true
    genotype.isHetNonRef shouldBe false
    genotype.getAlleles.map(_.getBaseString).toSeq should contain theSameElementsAs List("AAA", "AAAAAA")
    getStrGenotype(genotype) should contain theSameElementsInOrderAs Seq(1.0, 2.0)
  }

  it should "support the --min-cumulative-frequency option" in {
    val builder = new VariantContextSetBuilder(sampleNames=Seq.range(1, 7).map(i => s"S-$i"))
    val alleles = List("AAA", "AAAAAA", "AAAAAAAAA")

    // ref is most frequent (3), then allele #2 (2), then allele #1 (1)
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-1"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-2"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-3"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(2)), sampleName=Some("S-4"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(2)), sampleName=Some("S-5"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(1)), sampleName=Some("S-6"))

    // MCF 0.5 -> called homozygous ref
    {
      val outputs  = run(builder, minCumulativeFrequency=0.5)
      val variants = outputs.variants

      variants.length shouldBe 1
      variants.head.getGenotypes.length shouldBe 1
      val genotype = variants.head.getGenotype(0)
      genotype.isHom shouldBe true
      genotype.isHomRef shouldBe true
      genotype.getAlleles.map(_.getBaseString).toSeq should contain theSameElementsAs List(alleles(0), alleles(0))
      getStrGenotype(genotype) should contain theSameElementsInOrderAs Seq(1.0, 1.0)
    }

    // MCF 0.75 -> called heterozygous ref/alt
    {
      val outputs  = run(builder, minCumulativeFrequency=0.75)
      val variants = outputs.variants

      variants.length shouldBe 1
      variants.head.getGenotypes.length shouldBe 1
      val genotype = variants.head.getGenotype(0)
      genotype.isHet shouldBe true
      genotype.isHetNonRef shouldBe false
      genotype.getAlleles.map(_.getBaseString).toSeq should contain theSameElementsAs List(alleles(0), alleles(2))
      getStrGenotype(genotype) should contain theSameElementsInOrderAs Seq(1.0, 3.0)
    }

    // MCF 0.9 -> called triallelic
    {
      val outputs  = run(builder, minCumulativeFrequency=1.0)
      val variants = outputs.variants

      variants.length shouldBe 1
      variants.head.getGenotypes.length shouldBe 1
      val genotype = variants.head.getGenotype(0)
      genotype.isHet shouldBe true
      genotype.isHetNonRef shouldBe false
      genotype.getPloidy shouldBe 3
      genotype.getAlleles.map(_.getBaseString).toSeq should contain theSameElementsAs List(alleles(0), alleles(2), alleles(1))
      getStrGenotype(genotype) should contain theSameElementsInOrderAs Seq(1.0, 3.0, 2.0)
    }
  }

  it should "write no-calls when no genotypes overlap an STR" in {
    val builder  = new VariantContextSetBuilder(sampleNames=Seq.range(1, 7).map(i => s"S-$i"))
    val outputs  = run(builder, minCumulativeFrequency=0.0)
    val variants = outputs.variants
    variants.length shouldBe 1
    variants.head.getGenotypes.length shouldBe 1
    val genotype = variants.head.getGenotype(0)
    genotype.isNoCall shouldBe true
    genotype.hasAnyAttribute("STR_GT") shouldBe false
  }

  it should "support the --per-strand option" in {
    val builder = new VariantContextSetBuilder(sampleNames=Seq.range(1, 6).flatMap(i => Seq(s"S-$i/A", s"S-$i/B")))
    val alleles = List("AAA", "AAAAAA", "AAAAAAAAA")

    // ref is most frequent (3), then allele #2 (2)
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-1/A"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-1/B"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-2/A"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-2/B"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(0)), sampleName=Some("S-3/A")) // Missing 3/B!
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(2)), sampleName=Some("S-4/B")) // Missing 4/A!
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(2)), sampleName=Some("S-5/A"))
    builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List(alleles(2)), sampleName=Some("S-5/B"))

    val outputs = run(builder, minCumulativeFrequency=1.0, perStrand=true)
    val variants = outputs.variants

    variants.length shouldBe 1
    variants.head.getGenotypes.length shouldBe 1
    val genotype = variants.head.getGenotype(0)
    genotype.isHet shouldBe true
    genotype.isHetNonRef shouldBe false
    genotype.getAlleles.map(_.getBaseString).toSeq should contain theSameElementsAs List(alleles(0), alleles(2))
    getStrGenotype(genotype) should contain theSameElementsInOrderAs Seq(1.0, 3.0)
  }

  it should "write no-calls when there are no counts for the STR" in {
    val builder = new VariantContextSetBuilder(sampleNames=Seq.range(1, 7).map(i => s"S-$i"))
    val alleles = List("AAA", "AAAAAA", "AAAAAAAAA")

    // all no-calls
    Seq.range(1, 7).foreach { i =>
      builder.addVariant(refIdx=0, start=1, variantAlleles=alleles, genotypeAlleles=List.empty, sampleName=Some(s"S-$i"))
    }

    val outputs  = run(builder, minCumulativeFrequency=1.0)
    val variants = outputs.variants

    variants.length shouldBe 1
    variants.head.getGenotypes.length shouldBe 1
    val genotype = variants.head.getGenotype(0)
    genotype.isNoCall shouldBe true
    genotype.getAlleles.map(_.getBaseString).toSeq should contain theSameElementsAs List(Allele.NO_CALL_STRING, Allele.NO_CALL_STRING)
    getStrGenotype(genotype) should contain theSameElementsInOrderAs Seq(0.0, 0.0)
  }

  // TODO: add a test for aggregating genotype fields DP, GB, PDP, DSTUTTER, DFLANKINDEL from HipSTR
}
