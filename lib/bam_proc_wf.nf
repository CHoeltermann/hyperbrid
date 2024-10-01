#!/usr/bin/env nextflow

// sort, remove duplicates, and convert to bam

// code inspired from: https://unix.stackexchange.com/questions/512500/creating-a-list-of-files-in-bash-failed
process MERGE_BAM_by_ID {

	//container 'docker://biocontainers/bamtools:v2.5.1dfsg-3-deb_cv1'
	container 'docker://broadinstitute/gatk:4.6.0.0'

        label 'medium'

        input:
	// what do I take as input? a channel merged by SAMPLE_ID would probably be best!
        tuple val(SAMPLE_ID), val(RUN_ID), path(bam_files)

        output:
        tuple val(SAMPLE_ID), path("${SAMPLE_ID}_merged.bam")

	script:
	"""
	#bamtools merge -in ${bam_files.join(' -in ')} -out "${SAMPLE_ID}merged.bam"

	gatk MergeSamFiles \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		--INPUT ${bam_files.join(' --INPUT ')} \
		--OUTPUT "${SAMPLE_ID}_merged.bam"
	"""
}

process SORT_SAM {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

	input:
	tuple val(SAMPLE_ID), path(bam)

	output:
	tuple val(SAMPLE_ID), path("${SAMPLE_ID}_sort.bam"), path("*.bai")

	script:
	"""
	gatk SortSam \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-I $bam \
		-O "${SAMPLE_ID}_sort.bam" \
		--SORT_ORDER "coordinate"

	gatk BuildBamIndex \
		--java-options "-Xmx${task.memory.giga}G -XX:-UsePerfData" \
		-I "${SAMPLE_ID}_sort.bam"
	"""
}

process RMDUP {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

        input:
        tuple val(ID), path(bam), path(bai)

	publishDir { "${params.pubdir}/bam/" }, mode: 'copy'

	output:
	tuple(val(ID), path("${bam.baseName}_rmdup.bam"), path("${bam.baseName}_rmdup.bai"), emit: rmdup_bam)
	tuple(val(ID), path("${bam.baseName}_rmdup_metrics.txt"), emit: metrics)

        script:
        """
	gatk MarkDuplicates \
		--java-options "-Xmx${task.memory.giga-4}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-I "${bam.baseName}.bam" \
		-O "${bam.baseName}_rmdup.bam" \
		-M "${bam.baseName}_rmdup_metrics.txt" \
		--REMOVE_SEQUENCING_DUPLICATES true \
		--CREATE_INDEX true
	"""
}

process COV_CALC {

	label 'medium'

	container 'docker://broadinstitute/gatk:4.6.0.0'

	input:
	path ref_fa
	tuple val(ID), path(bam), path(bai)

	publishDir { "${params.pubdir}/cov/" }, pattern: '*_picard_cov.txt', mode: 'copy'

	output:
	path "${ID}_picard_cov.txt", emit: qc
	tuple(val(ID), env(COV), path(bam), path(bai), emit: bam)

	shell:
	'''
	gatk CollectWgsMetrics \
		--java-options "-Xmx!{task.memory.giga}G -XX:-UsePerfData" \
		--TMP_DIR ./tmp \
		-I !{bam} \
		-O "!{ID}_picard_cov.txt" \
		-R !{ref_fa} \
		--VALIDATION_STRINGENCY SILENT \
		--READ_LENGTH 500 \
		--USE_FAST_ALGORITHM true

	# extract mean cov as value here
	COV=$( sed -n -e 8p "!{ID}_picard_cov.txt" | cut -f 2 )
	'''
}

workflow BamProc {

	take:
	ref_fa
	ref_fai
	ref_dict
	reads
	bam_file

	main:
	merged_bam_ch = bam_file.groupTuple(by: 0)
	merged_bam_ch.view()
	MERGE_BAM_by_ID(merged_bam_ch)
	
	SORT_SAM(MERGE_BAM_by_ID.out)

	RMDUP(SORT_SAM.out)
	COV_CALC(ref_fa, RMDUP.out.rmdup_bam)

	emit:
	RMDUP_BAM = COV_CALC.out.bam
	RMDUP_PICARD_METRICS = RMDUP.out.metrics
	COV = COV_CALC.out.qc
}

