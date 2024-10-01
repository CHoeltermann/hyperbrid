#!/usr/bin/env nextflow

// Do GWAS GQ & analysis

// 7 steps

// 1. Missingness of SNPs and individuals
// 2. sex discrepancy
// 3. minor allele frequency
// 4. deviations from hardy weinberg equillibrium
// 5. heterozygosity rate
// 6. relatedness
// 7. ethnic outliers -- dont know if this is applicable to my monkeys


process {

	label ''

	container 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'

	publishDir { "${params.pubdir}/" }, mode: 'copy'

	input:

	output:

	script:
	"""
	"""

}

process {

	label ''

	container 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'

	publishDir { "${params.pubdir}/" }, mode: 'copy'

	input:

	output:

	script:
	"""
	"""

}

process {

	label ''

	container 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'

	publishDir { "${params.pubdir}/" }, mode: 'copy'

	input:

	output:

	script:
	"""
	"""

}

workflow {



}

