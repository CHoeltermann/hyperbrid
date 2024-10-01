#!/usr/bin/env nextflow

// following snpEFF v4.3 (Cingolani et al. 2012) best-practice pipeline for annotation
// https://pcingola.github.io/SnpEff/snpeff/inputoutput/#effect-prediction-details
// "SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes)"

// I need an option to decide wether to build a new SNP anno reference, or wether to search in the snpEff for an existing one.
// How do I build one?

// accessing existing ones:
/*
-- from inside the singularity container, it's weirder, bc the path is needed
singularity shell https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2 --no-home
cd /usr/local/share/snpeff-5.1-2/
java -jar snpEff.jar databases

*/

// if  I aligned against Mmul10, I can download it:
process SNP_DB {

	container 'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2'

	input:
	path ref_fa
	path ref_gtf
	path ref_cdna_fa

	output:
	path "./*/", emit: snpEff_db

	script:
	if ( params.snpEff_download == true )
	"""
	java \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-c /usr/local/share/snpeff-5.1-2/snpEff.config \
		-noLog \
		-t \
		-noShiftHgvs \
		-jar /usr/local/share/snpeff-5.1-2/snpEff.jar \
		download $params.snp_ref
	"""

	else
	"""
	# so, Apparently the snpEff.config from the container needs to be updated here. I guess thats possible, but kinda extensive.
	# a new entry must be added there
	#
	# example (add this to snpEff.config, if i had a new mouse genome):
	# # Mouse genome, version mm37.61
	# mm37.61.genome : Mouse
	#
	# supply gtf.gz, fa.gz as input pathjs for this channel. This should link the existing files; and the "build" command should recognize them as "present" in the wd.
	# give genome name from above (config file) as an option and it should work.
	
	# protein CDS are used by snpEff for sanity check / consistency check
	
	#java \
	#	-Xmx"${process.memory}" \
	#	-c /usr/local/share/snpeff-5.1-2/snpEff.config \
	#	-noLog \
	#	-t \
	#	-noShiftHgvs \
	#	-jar /usr/local/share/snpeff-5.1-2/snpEff.jar \
	#	build mm37.61
	"""
}

process SNPeffANNO {

	container 'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2'

        publishDir { "${params.pubdir}/vcf/" }, mode: 'copy'

	input:
	path db
	path vcf

	output:
	tuple val(ID), path("${vcf.baseName}_anno.vcf")

	script:
	"""
	java \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-c /usr/local/share/snpeff-5.1-2/snpEff.config \
		-noLog \
		-t \
		-noShiftHgvs \
		-jar /usr/local/share/snpeff-5.1-2/snpEff.jar \
		anno > "${vcf.baseName}_anno.vcf" 
	"""
}

workflow SNPanno {

	take:
	ref_fa
	ref_gtf
	ref_cdna_fa
	vcf

	main:
	SNP_DB(ref_fa, ref_gtf, ref_cdna_fa)
	db_ch = SNP_DB.out.snpEff_db.collect()
	SNPeffANNO(db_ch, vcf)

	emit:
	SNPeffANNO.out.anno_vcf
}

