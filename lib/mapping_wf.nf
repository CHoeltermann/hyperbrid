#!/usr/bin/env nextflow
//https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently

process IDX_RG {

        label 'medium'

        container 'https://depot.galaxyproject.org/singularity/bwa:0.7.16--pl5.22.0_0'

        input:
        path ref_fa

        output:
        path "*", emit: idx 
        path ref_fa, emit: ref_fa

        script:
        """
        bwa index ${ref_fa}
        """
}

// read cleaning

process FASTP {

        label 'small'

	container 'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h5f740d0_3'

        input:
        tuple val(SAMPLE_ID), val(RUN_ID), path(R1), path(R2)

        output:
        tuple(val(SAMPLE_ID), val(RUN_ID), path("${RUN_ID}_clean_1.fastq.gz"), path("${RUN_ID}_clean_2.fastq.gz"), emit: reads)
        path("${RUN_ID}_fastp.json"), emit: info

        script:
        """
        fastp \
                -i $R1 \
                -I $R2 \
                -o "${RUN_ID}_clean_1.fastq.gz" \
                -O "${RUN_ID}_clean_2.fastq.gz" \
                -j "${RUN_ID}_fastp.json"
        """
}

// BAM requirements by GATK/ Picard: https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats
process ALIGN {

        label 'intense'

        container 'docker://michaelfranklin/bwasamtools:0.7.17-1.10'
        // container 'https://depot.galaxyproject.org/singularity/bwa:0.7.16--pl5.22.0_0'

        input:
        tuple val(SAMPLE_ID), val(RUN_ID), path(R1), path(R2)
        path ref_files_idx
        path ref_fa

        output:
        tuple(val(SAMPLE_ID), val(RUN_ID), path("${RUN_ID}.bam"), emit: bam)

        script:
        """
        bwa mem \
                -t ${task.cpus} \
                -M \
		-R "@RG\\tID:${RUN_ID}\\tSM:${SAMPLE_ID}\\tLB:${SAMPLE_ID}-${RUN_ID}\\tPL:${params.seq_platform}" \
                $ref_fa \
                $R1 \
                $R2 \
                | samtools view \
                -b \
		-h \
                -T $ref_fa > "${RUN_ID}.bam"
        """
}

workflow Mapping {

	take:
	ref_fa
	ch_reads

	main:
	IDX_RG(ref_fa)
	// make value ch
	idx_ch = IDX_RG.out.idx.collect()

	FASTP(ch_reads)
	ALIGN(FASTP.out.reads, idx_ch, ref_fa)

	emit:
	bam_aln = ALIGN.out.bam
	clean_fq = FASTP.out.reads
	fastp_qc = FASTP.out.info
}

