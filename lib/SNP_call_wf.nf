#!/usr/bin/env nextflow

// Nice tutorial:
// https://hpc.nih.gov/training/gatk_tutorial/workflow-overview.html
// https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/

// I can use Picard from GATK; it is shipped within GATK in the newer versions. This just means I can get away with one less container, i guess
// remade after best practices workflow: https://github.com/gatk-workflows/gatk4-germline-snps-indels

// 2. MergeGVCFs

// --- SNP calling
// realignment of indel polymorphisms with RealignerTargetCreator & IndelRealigner
// recalibrate QC scores, correct previously called variants (HaplotypeCaller)
// 

// --- Genotype calling
// GVCF for each individual using HaplotypeCaller
// GenotypeGVCFs based method with the “includeNonVariantSites” flag, to get the population VCF file, including all confident sites
// “SelectVariants” to exclude indels
// split data into variant & nonvariant sites
//  “VariantFiltration” was applied to exclude potential false positive variant calls with the following criteria: “filterExpression QD < 2.0 || FS > 60.0 || MQ < 40.0 || ReadPosRankSum < --8.0 || MQRankSum < 12.5” and “genotypeFilterExpression DP < 4.0”
// sites were removed if there was an “N” in the reference sequence or the site spanned an indel plus a buffer of 3 bp in both directions and the site included >10% missing genotypes
// To obtain the genotype file for subsequent analyses, a PERL script was used to transfer the VCF format to genotype format (e.g. AA, AT) and degenerate bases format (e.g. ‘M’ = ‘AC’) for all langurs and then again to generate final genotypes in VCF format.

// --- Annotation of SNPs
// snpEFF v4.3 best practise pipeline - with custom Tfra_2.0 dataset (gtf/gff3 from NCBI, supp. figure 15, supp tables 1 and 24)


// https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode
// HaploCtypeCaller implicitly expects a .fai file! I'll need to forward it as input, as well.

// this is highly recommended by gatk, but needs a database of SNPs, I think
// did Liye do this? Should I do this?
// In best practises of GATK, its labelled as "optional". 

process GATK_QC {

        label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

        input:
        path ref_fa
        path ref_fai
        path ref_dict
        tuple val(ID), path(bam), path(bai)

        output:
	path "${bam.baseName}_sumMetricsGATK.txt", emit: metrics
	// tuple(path("${bam.baseName}_sumMetricsGATK.txt"), path(path "${bam.baseName}_insertMetricsGATK.txt"), path("${bam.baseName}_sumMetricsGATKhist.pdf"), emit: metrics)

        script:
        """
	# do this after dup removal, or before?
        gatk CollectAlignmentSummaryMetrics \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
                -R $ref_fa \
                -I $bam \
                -O "${bam.baseName}_sumMetricsGATK.txt" # output; .txt file
        
        #gatk CollectInsertSizeMetrics \
	#	--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
        #        -I $bam \
        #        -O "${bam.baseName}_insertMetricsGATK.txt" \
        #        -H "${bam.baseName}_sumMetricsGATKhist.pdf" \
        #        -M 0.05
        """
}

process HAPLOCALL {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

	input:
	path ref_fa
	path ref_fai
	path ref_dict
	tuple val(ID), path(bam), path(bai)

	output:
	tuple(val(ID), path("${bam.baseName}_raw.g.vcf.gz"), path("${bam.baseName}_raw.g.vcf.gz.tbi"), emit: vcf_raw)
        tuple(val(ID), path(bam), path(bai), emit: bam_bai)

	script:
	"""
	# run per-sample; HC should not be called on more than 100 samples - parallelizing solves this

	gatk HaplotypeCaller \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		--reference $ref_fa \
		--input $bam \
		-output "${bam.baseName}_raw.g.vcf.gz" \
		-ERC GVCF 
	"""
}

/*
// since this step leverages population-wide information, it must be done on a per-population basis & cant be parallelized
process MERGEGVCF {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

	input:
	path ref_fa
	path ref_fai
	path ref_dict
	tuple val(ID), path(vcf), path(vcf_tbi)

	output:
	tuple(val(ID), path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi"), emit: gvcf)

	script:
	"""
	# I think the .collect() does not work like this ... it only finds one file.
	#find . -name "*.g.vcf.gz" > gvcf.list
	ls $vcf > gvcf.list

	gatk CombineGVCFs \
		--TMP_DIR ./tmp \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		-R $ref_fa \
		--variant "gvcf.list" \
		-O "cohort.g.vcf.gz"
	"""
}
*/

process GENOGVCF {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

	input:
	path ref_fa
	path ref_fai
	path ref_dict
	tuple val(ID), path(gvcf), path(vcf_idx)

	output:
	tuple(val(ID), path("genotyped.vcf.gz"), path("genotyped.vcf.gz.tbi"), emit: vcf_geno)

	script:
	"""
	gatk GenotypeGVCFs \
		--TMP_DIR ./tmp \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		-R $ref_fa \
		-V $gvcf \
		-O "genotyped.vcf.gz" \
		--include-non-variant-sites
	"""
}

// I'll likely need to change the parameters often so I should make it accessible
process VFILTER {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'
	
	input:
	path ref_fa
	path ref_fai
	path ref_dict
	tuple val(ID), path(vcf_geno), path(vcf_geno_tbi)

	output:
	tuple(val(ID), path("${vcf_geno.baseName}_filtered.vcf.gz"), path("${vcf_geno.baseName}_filtered.vcf.gz.tbi"), emit: vcf_filtered)

	script:
	"""
	gatk VariantFiltration \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-R $ref_fa \
		-V $vcf_geno \
		-O "${vcf_geno.baseName}_filtered.vcf.gz" \
		--verbosity ERROR \
		-filter ${params.filter.expr1} --filter-name ${params.filter.name1} \
		-filter ${params.filter.expr2} --filter-name ${params.filter.name2} \
		-filter ${params.filter.expr3} --filter-name ${params.filter.name3} \
		-filter ${params.filter.expr4} --filter-name ${params.filter.name4} \
		-filter ${params.filter.expr5} --filter-name ${params.filter.name5} \
		-filter ${params.filter.expr6} --filter-name ${params.filter.name6}

	# At this point I also need to filter for the following:
	# If there is "N" in the reference sequence
	# "or the site spanned an indel plus a buffer of 3 bp in both directions and the site included >10% missing genotypes"
	# How do I do this most simply? By using gatk again, can it do that? BCFtools? Filter with bash/python/R script myself? 
	"""
}

// extract SNPs/ Indels
// exclude filtered variants (file has FAIL/PASS variable, generated in VFILTER)
// split also into variant & nonvariant sites
process SELECT {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'
	
	publishDir { "${params.pubdir}/vcf/" }, mode: 'copy'

	input:
	path ref_fa
	path ref_fai
	path ref_dict
	tuple val(ID), path(vcf_filtered), path(vcf_filtered_idx)

	output:
	tuple(val(ID), path("${vcf_filtered.baseName}_snp.vcf.gz"), path("${vcf_filtered.baseName}_snp.vcf.gz.tbi"), emit: snp)
	tuple(val(ID), path("${vcf_filtered.baseName}_indel.vcf.gz"), path("${vcf_filtered.baseName}_snp.vcf.gz.tbi"), emit: indel)

	script:
	"""
	gatk SelectVariants \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-R $ref_fa \
		-V $vcf_filtered \
		--select-type-to-include SNP \
		--exclude-filtered \
		--create-output-variant-index true \
		-O "${vcf_filtered.baseName}_snp.vcf"

	# only need SNPs, can make this step optional with default to "false"
	gatk SelectVariants \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-R $ref_fa \
		-V $vcf_filtered \
		--select-type-to-include INDEL \
		--exclude-filtered \
		--create-output-variant-index true \
		-O "${vcf_filtered.baseName}_indel.vcf"

	# I'll also need to filter out those that failed the genotype filters
	# grep can be used, but maybe SelectVariants already has an option for that?
	"""
}

workflow SNPcall {

        take:
	ref_fa
	ref_fai
	ref_dict
	bam_bai

        main:
	HAPLOCALL(ref_fa, ref_fai, ref_dict, bam_bai)
	GATK_QC(ref_fa, ref_fai, ref_dict, HAPLOCALL.out.bam_bai)
	
	//vcf_paths_ch = HAPLOCALL.out.vcf_raw.collect()

	GENOGVCF(ref_fa, ref_fai, ref_dict, HAPLOCALL.out.vcf_raw)
	VFILTER(ref_fa, ref_fai, ref_dict, GENOGVCF.out.vcf_geno)
	SELECT(ref_fa, ref_fai, ref_dict, VFILTER.out.vcf_filtered)

        emit:
	QC = GATK_QC.out.metrics
	SNP = SELECT.out.snp
	INDEL = SELECT.out.indel
}

