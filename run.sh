#!/bin/bash
#NXF_VER=24.04.2

nextflow run ./main.nf \
	-profile standard \
	-resume \
	--samplesheet "samplesheet.csv" \
	--chr_map "chr_map_reduced.tsv" \
	--refg_fa "reference_genome.fna.gz" \
        --refg_gtf "reference.gtf.gz" \
        --snp_ref "Mmul_10.105" \
	--max_cpus 8 \
	--max_memory '24.GB' #\


