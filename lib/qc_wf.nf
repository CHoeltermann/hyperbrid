#!/usr/bin/env nextflow

process FASTQC {

        label 'medium'

        container 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'

        publishDir { "${params.pubdir}/fastqc/" }, pattern: '*.{zip,html}',  mode: 'copy'

        input:
        tuple val(SAMPLE_ID), val(RUN_ID), path(R1), path(R2)
        tuple val(SAMPLE_ID), val(RUN_ID), path(bam)

        output:
        path "*.{zip,html}"

        script:
        """
        fastqc $R1 $R2 -q -o .
        fastqc $bam -q -o .
        """
}

// https://github.com/scchess/Qualimap
process QUALIMAP {

        label 'medium'

        container 'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--hdfd78af_2'

        input:
        tuple val(ID), path(bam), path(bai)

        publishDir { "${params.pubdir}/qualimap/" }, mode: 'copy'

        output:
        path "*"

        script:
        """
        # bam needs to be sorted by chromosome coordinates!
        qualimap bamqc -bam $bam #\
                #--java-mem-size="${task.memory.giga}G"
        """ 
}

// MultiQC does not really work as expected; needs work.
workflow QualiWorker {

	take:
		ch_reads
		bam_ch

	main:
		FASTQC(ch_reads)
		QUALIMAP(bam_ch)

	emit:
		fastqc = FASTQC.out
		qualimap = Qualimap.out
}

