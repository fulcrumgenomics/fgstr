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

import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.FgBioDef.IteratorToJavaCollectionsAdapter
import com.fulcrumgenomics.cmdline.ClpGroups
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.str.cmdline.FgStrTool
import com.fulcrumgenomics.str.vcf.StrInterval._
import com.fulcrumgenomics.util.{ProgressLogger, Rscript}
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.ReferenceSequenceFileWalker
import htsjdk.samtools.util.{Interval, IntervalList, StringUtil}
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.collection.mutable.ListBuffer
import scala.util.Failure

@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Calls a diploid genotype for individually called duplex molecules.
    |
    |The input VCF should be produce by joint calling each duplex source molecule as though they are separate samples.
    |This can be achieved by transforming the output of fgbio's `GroupReadsByUmi` with `ReadGroupPerDuplexMolecularId`,
    |and then running it through HipSTR, to be produce a multi-sample VCF.  The following options should be used when
    |running HipSTR:
    |
    |```
    |  --regions <path-to-regions>.bed \
    |  --output-gls \
    |  --output-pls \
    |  --min-reads 1 \
    |  --use-unpaired \
    |  --no-rmdup \
    |  --hap-chr-file <path-to-a-list-of-all-chromosomes>.txt \
    |  --max-str-len 150 \
    |  --max-flank-haps 10 \
    |  --filter-flank-haps
    |```
    |
    |NB: `--max-flank-haps` and `--filter-flank-haps` require the following branch:
    |  https://github.com/tfwillems/HipSTR/pull/30
    |
    |An interval list of specifying the same regions as input to HipSTR should be given.  The name field should contain
    |a comma list of values as follows:
    |  1. the repeat unit length (ex. `3` for the tri-nucleotide repeat `TCATCATCATCA`).
    |  2. the number of repeat units (ex. `4` for the tri-nucleotide repeat `TCATCATCATCA`).
    |  3. the name of the STR (ex. D1S1656)
    |  4. optionally, the expected (known or truth) number of repeat units for allele #1
    |  5. optionally, the expected (known or truth) number of repeat units for allele #2
    |
    |The output will contain a single sample with a diploid genotype call.  The INFO field will have list the allele
    |frequencies for all alleles observed, not just those in the final genotype.
    |
    |A PDF with plots will be output, containing a plot per STR showing the counts observed repeat unit length.  If the
    |expected number of repeat units are given, a horizontal dotted line will be plotted for each expected allele.
    |
    |If the `--per-strand` option is set, then the input should contain a sample per duplex source molecule and strand
    |(i.e. A and B strand respectively).  The calls from both strands will be compared for perfect agreement.  If only
    |one strand is present, then that call is retained.  The produced per duplex source molecule genotype calls are then
    |treated normally.
  """)
class StrGenotypeDuplexMolecules
(
  @arg(flag='i', doc="Input VCF of variant calls from HipSTR.") val input: PathToVcf,
  @arg(flag='o', doc="Prefix for all output files.") val output: PathPrefix,
  @arg(flag='r', doc="Reference sequence fasta file.") val ref: PathToFasta,
  @arg(flag='l', doc="Interval list with the STR regions.") val intervals: PathToIntervals,
  @arg(flag='S', doc="The output sample name") val sampleName: String = "Sample",
  @arg(flag='s', doc="Expect a sample per-strand of a duplex molecule") val perStrand: Boolean = false,
  @arg(flag='m', doc="Minimum minor allele frequency to call a heterozygous genotype") val minimumMaf: Double = 0.2,
  private val skipPlots: Boolean = false // for not plotting in tests
) extends FgStrTool with LazyLogging {

  private val ScriptPath = "com/fulcrumgenomics/str/vcf/StrGenotypeDuplexMolecules.R"

  Io.assertReadable(input)
  Io.assertReadable(intervals)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    def f(ext: String): FilePath = PathUtil.pathTo(output + ext)

    val vcfIn        = new VCFFileReader(input.toFile, false)
    val strIntervals = StrInterval.loadIntervals(this.intervals)
    val ref          = new ReferenceSequenceFileWalker(this.ref.toFile)
    val dict         = ref.getSequenceDictionary
    val vcfOut       = toVcfWriter(header=vcfIn.getFileHeader, sampleName=sampleName, dict=dict, output=f(".vcf.gz"))
    val infoWriter   = Io.toWriter(f(".info.txt"))
    val progress     = ProgressLogger(logger=logger, noun="variants")


    verifySampleNames(vcfIn.getFileHeader)

    strIntervals.foreach { str =>
      var numOverlapping = 0
      // go through all the variant calls that overlap
      vcfIn
        .query(str.chrom, str.start, str.end)
        .map(fromSingleStrand)
        .foreach { ctx =>
          numOverlapping += 1

          val ctxBuilder = new VariantContextBuilder(ctx)
          val genotypes  = ctx.getGenotypes
          val alleles    = ctx.getAlleles.toSeq

          // get the count per allele
          val counts: Seq[Int] =  alleles.map { allele => genotypes.map { genotype => genotype.countAllele(allele) }.sum }
          val totalCounts: Int = counts.sum

          // get the new genotype
          val allCalls = StrAllele.toCalls(str, ctx, counts)

          // naively choose the two most frequent genotypes
          val genotypeCalls = {
            allCalls.take(2) match {
              case Seq()             =>
                val noCall = StrAllele(Allele.NO_CALL, 0, 0)
                Seq(noCall, noCall)
              case Seq(call)         =>
                Seq(call, call)
              case Seq(call1, call2) =>
                //val maf: Double = call2.count / (call1.count + call2.count).toDouble
                val maf: Double = call2.count / totalCounts.toDouble
                if (maf < minimumMaf) Seq(call1, call1) else Seq(call1, call2)
            }
          }
          val genotypeBuilder = new GenotypeBuilder(sampleName, genotypeCalls.map(_.allele).toIterator.toJavaList)

          addFormatFields(
            ctx        = ctx,
            builder    = genotypeBuilder,
            calls      = genotypeCalls,
            unitLength = str.unitLength
          )

          // build the variant context
          ctxBuilder.genotypes(genotypeBuilder.make())

          // write the outputs
          vcfOut.add(ctxBuilder.make())
          infoWriter.write(str.toLongString(allCalls) + "\n")
          progress.record(ctx.getContig, ctx.getStart)
        }
      // if not calls overlapped, then output a missing genotype
      if (numOverlapping == 0) {
        val referenceSequence = ref.get(dict.getSequenceIndex(str.chrom))
        val referenceBases    = referenceSequence.getBases
        val alleleBases       = StringUtil.bytesToString(referenceBases, str.start-1, str.length)
        val allele            = Allele.create(alleleBases, true)
        val ctxBuilder        = new VariantContextBuilder(this.getClass.getSimpleName, str.chrom, str.start, str.end, Iterator(allele).toJavaList)
        val genotypeBuilder   = new GenotypeBuilder(sampleName, Iterator(Allele.NO_CALL, Allele.NO_CALL).toJavaList)
        val ctx               = ctxBuilder.genotypes(genotypeBuilder.make()).make()
        // write the outputs
        vcfOut.add(ctxBuilder.make())
        infoWriter.write(str.toLongString(Seq(StrAllele.NoCall, StrAllele.NoCall)) + "\n")
        progress.record(ctx.getContig, ctx.getStart)
      }
    }

    vcfIn.safelyClose()
    vcfOut.close()
    infoWriter.close()

    if (!skipPlots && progress.getCount > 0) {
      Rscript.execIfAvailable(ScriptPath, f(".info.txt").toString, f(".pdf").toString) match {
        case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
        case _ => Unit
      }
    }
    else {
      logger.warning("No variants outputted, skipping plots")
    }
  }

  private def fromSingleStrand(ctx: VariantContext): VariantContext = {
    if (!this.perStrand) return ctx

    val ctxBuilder = new VariantContextBuilder(ctx)

    val pairs = ctx.getGenotypes.toSeq.groupBy { genotype =>
      val miStrand = {
        val name  = genotype.getSampleName
        val index = name.lastIndexOf('-')
        require(index > 0, s"Sample name doesn't look like a duplex id: $name")
        name.substring(index + 1)
      }
      val index = miStrand.lastIndexOf('/')
      require(index > 0, s"Molecular id doesn't look like a single strand id: $miStrand")
      miStrand.substring(0, index)
    }

    val genotypes = pairs.flatMap { case (mi, gs) =>
      require(gs.length == 2 || gs.length == 1, s"Expected one or two genotypes (per-strand), found ${gs.length}")
      if (gs.map(_.getGenotypeString).toSet.size == 1) { // they should all agree
        val name  = gs.head.getSampleName
        val index = name.lastIndexOf('/')
        val genotypeBuilder = new GenotypeBuilder(gs.head) // TODO: merge format fields, for now keep the first
        genotypeBuilder.name(name.substring(0, index))
        Some(genotypeBuilder.make())
      }
      else {
        None
      }
    }

    ctxBuilder.genotypes(genotypes.toIterator.toJavaList)
    ctxBuilder.make()
  }

  /** Adds various format fields that can be easily aggregated across genotypes (ex. DP, GB, PDP, DSTUTTER, DFLANKINDEL) */
  private def addFormatFields(ctx: VariantContext, builder: GenotypeBuilder, calls: Seq[StrAllele], unitLength: Int): Unit = {
    val alleles = calls.map(_.allele)

    var dp: Int = 0
    val gbs: ListBuffer[Option[Int]] = ListBuffer(None, None)
    val pdps: ListBuffer[Double] = ListBuffer(0.0, 0.0)
    var dStutter: Int = 0
    var dFlankIndel: Int = 0

    ctx.getGenotypes.foreach { genotype =>
      val genotypeAlleles = genotype.getAlleles
      require(genotypeAlleles.length == 1, "Found more than one allele, expected a haploid call")
      val genotypeAllele = genotypeAlleles.get(0)

      alleles.zipWithIndex.find(_._1 == genotypeAllele).foreach { case (allele, alleleIndex) =>

        // DP
        dp += genotype.getAttributeAsInt("DP", 0)

        // GB
        if (genotype.hasAnyAttribute("GB")) {
          val gbAttribute: Int = genotype.getAttributeAsInt("GB", 0)
          gbs(alleleIndex) match {
            case None    => gbs(alleleIndex) = Some(gbAttribute)
            case Some(v) => require(v == gbAttribute, s"Mismatching GB field, expected '$v', found '$gbAttribute', for allele $genotypeAllele")
          }
        }

        // PDP
        val attributeIndex = genotypeAlleles.indexOf(allele)
        if (genotype.hasAnyAttribute("PDP")) {
          val pdpAttibutes: List[Double] = genotype.getAnyAttribute("PDP").asInstanceOf[java.util.List[java.lang.String]].toList.map(_.toDouble)
          pdps(alleleIndex) += pdpAttibutes(attributeIndex)
        }

        // DSTUTTER
        if (genotype.hasAnyAttribute("DSTUTTER")) {
          dStutter += genotype.getAttributeAsInt("DSTUTTER", 0)
        }

        // DFLANKINDEL
        if (genotype.hasAnyAttribute("DFLANKINDEL")) {
          dFlankIndel += genotype.getAttributeAsInt("DFLANKINDEL", 0)
        }
      }
    }

    // set them
    builder.attribute("DP", dp)
    if (gbs.forall(_.isDefined)) builder.attribute("GB", gbs.flatten.toIterator.toJavaList)
    builder.attribute("PDP", pdps.toIterator.toJavaList)
    builder.attribute("DSTUTTER", dStutter)
    builder.attribute("DFLANKINDEL", dFlankIndel)
    builder.attribute("STR_GT", calls.map(_.repeatLength / unitLength.toDouble).toIterator.toJavaList)
  }

  /** Creates the output VCF writer */
  private def toVcfWriter(header: VCFHeader, sampleName: String, dict: SAMSequenceDictionary, output: PathToVcf): VariantContextWriter = {
    // Developer notes:
    // - at a some point in the future, it would be nice to remove all of the INFO and FORMAT definitions we no longer
    //   use
    val headerLines: util.Set[VCFHeaderLine] = new util.HashSet[VCFHeaderLine](header.getMetaDataInSortedOrder)
    headerLines.add(new VCFFormatHeaderLine("STR_GT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "The STR genotype (number of units)"))
    val outputHeader = new VCFHeader(headerLines, Iterator(sampleName).toJavaSet)
    if (outputHeader.getSequenceDictionary == null) outputHeader.setSequenceDictionary(dict)
    val builder = new VariantContextWriterBuilder()
      .setOutputFile(output.toFile)
      .setReferenceDictionary(header.getSequenceDictionary)
      .setOption(Options.INDEX_ON_THE_FLY)
      .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, true)
    val writer = builder.build
    writer.writeHeader(outputHeader)
    writer
  }

  /** Verifies that all the input sample names to be the same after removing the trailing "mid-<mid>". */
  private def verifySampleNames(header: VCFHeader): String = {
    val sampleNames = header
      .getSampleNamesInOrder
      .map { name =>
        val index = name.lastIndexOf('-')
        require(index > 0, s"Sample name '$name' does not have a trailing '-<mid>'")
        name.substring(0, index)
      }.toSet.toList
    sampleNames match {
      case Nil         => fail("No sample names found!")
      case name :: Nil => name
      case names       => fail(s"More than one sample name found: ${names.mkString(", ")}")
    }
  }
}