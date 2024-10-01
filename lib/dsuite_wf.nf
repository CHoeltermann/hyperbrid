#!/usr/bin/env nextflow

// use Dsuite to determine if further analysis makes sense; and if yes, for which samples!

process DSUITE {

        container 'docker://zjnolen/dsuite:0.5-r52'

        input:
        tuple val(ID), path(vcf), path(vcf_idx), path(sets)

        publishDir { "${params.pubdir}/dsuite" }, mode: 'copy'

        // use Dval to filter channel to Dinvestigate; maybe even with a params.Dval from nextflow.config
        output:
        tuple val(ID), path("BBAA.txt"), path("Dmin.txt"), val(Dval)

        script:
        """ 
        Dsuite Dtrios \
                --JKnum 20 \
                $vcf \
                $sets
        """ 
}

process DSUITE_TREE {

        container 'docker://zjnolen/dsuite:0.5-r52'

        input:
        tuple val(ID), path(vcf), path(vcf_idx), path(sets)
        path tree

        publishDir { "${params.pubdir}/dsuite" }, mode: 'copy'

        output:
        tuple val(ID), path("*BBAA.txt"), path("*Dmin.txt"), path("*tree.txt")

        script:
        """ 
        Dsuite Dtrios \
                --tree=$tree \
                --JKnum 20 \
                --ABBAclustering \
                $vcf \
                $sets
        """
}


process D_INV {

        container 'docker://zjnolen/dsuite:0.5-r52'

        input:
        tuple val(ID), path(vcf), path(vcf_idx), path(sets)

        publishDir { "${params.pubdir}/dsuite" }, mode: 'copy'

        output:

        script:
        """
	# make test_trio.txt
	# use bash; I can get trios that are interesting - how???
        Dsuite Dinvestigate \
                [OPTIONS] \
                $vcf \
                $sets \
                test_trios.txt
        """
}

// DSUITE -- https://github.com/millanek/Dsuite
workflow Dsuite {

        take:
        vcf

        main:
        if ( params.treefile_exists ) {
                tree_ch = Channel.fromPath(params.treefile)
                DSUITE_TREE(vcf, tree_ch)
        }
        else {
                DSUITE(vcf)
        }
        // filter for D statistic, code logic to decide which to Dinvestigate

        emit:
        DSUITE.out
}

