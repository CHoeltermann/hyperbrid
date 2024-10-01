#!/usr/bin/env nextflow
// ADMIXTURE

process DSUITE {

	// needs a population/species map (a file)
	// get from reads file ...
        container 'docker://zjnolen/dsuite:0.5-r52'

        input:
	tuple val(ID), path(vcf), path(vcf_idx)
	path sets

        publishDir { "" }, mode: 'copy'

        output:
	tuple val(ID), path("BBAA.txt"), path("Dmin.txt"), path("tree.txt")

        script:
	// params.groups
	// newick tree files are optional
        """
	Dsuite Dtrios \
		$vcf \
		$sets
        """
}

workflow ADMIX {

        take:
	vcf
	sets

        main:
	DSUITE(vcf, sets)

        emit:
	DSUITE.out
}

