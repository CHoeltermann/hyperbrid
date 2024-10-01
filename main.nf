#!/usr/bin/env nextflow

// Params
// these can be overwritten on the cmd line, like: --reads ...

log.info """\
    ===================================
     H Y P E R - N F   P I P E L I N E
    ===================================
	samplesheet       : ${params.samplesheet}
        output to         : ${params.pubdir}
        started at        : ${workflow.start}
        config files      : ${workflow.configFiles}

	--
	minimum coverage  : ${params.min_cov}
	downsampling      : ${params.downsampling}
	outgroup (tree)   : ${params.outgroup}
         """
         .stripIndent()

// import from modules

include { PCA         } from './lib/pca_wf.nf'
include { ADMIX       } from './lib/admixture_wf.nf'
include { Mapping     } from './lib/mapping_wf.nf'
include { BamProc     } from './lib/bam_proc_wf.nf'
include { SNPcall     } from './lib/SNP_call_wf.nf'
include { SNPanno     } from './lib/SNP_anno_wf.nf'
include { QualiWorker } from './lib/qc_wf.nf'

// Processes

// some prep steps:

// test input files to see of FASTQ format is ok

process FAIDX {

	label 'medium'

	container 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0'

	input:
	path ref_fa

	output:
	path "${ref_fa.baseName}", emit: ref_fa
	path "*.fai", emit: ref_fai

	script:
	"""
	gzip -d -k $ref_fa
	samtools faidx "${ref_fa.baseName}" 
	"""
}

process DICT {

	label 'medium'

        container 'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.0--hdfd78af_0'

        input:
        path ref_fa

        output:
        path "*.dict", emit: ref_dict

        script:
        """
        gatk CreateSequenceDictionary \
                -R $ref_fa
        """
}

process DOWNSAMPLING {

	label 'medium'

	container 'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.0--hdfd78af_0'

	input:
	tuple val(ID), val(cov_val), path(bam), path(bai)

	output:
	path "*_downsampled.bam"

	// I need to calculate P based on the coverage I have, and the coverage I want.
	shell:
	'''
	PER=params.min_cov/!{cov_val}
	
	gatk DownsampleSam \
		-I !{bam} \
		-O !{bam.baseName}_downsampled.bam \
		-P $PER
	'''
}

// take here all vcf files - SNP filtered and everything
process VCF_MERGE {

	label 'medium'

	publishDir { "${params.pubdir}/vcf/" }, mode: 'copy'

        container 'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_1'

        input:
        path(snp)

	output:
	tuple path("merged.vcf.gz"), path("*.vcf.gz.tbi")

	// merge and index?!
	script:
	"""
	bcftools merge \
		--force-samples \
		--gvcf \
		--output-type z \
		--write-index="tbi" \
		--output "merged.g.vcf.gz" \
		*vcf.gz

#	bcftools view "merged.g.vcf.gz" \
#		--regions ^X,^Y
#		--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21 \
#		--output-type z \
#		--output "${snp.baseName}_autosomal.vcf.gz"
	"""
}

process RENAME_CHR_VCF {

	input:
	tuple path(vcf), path(vcf_idx)

	output:
	tuple path("merged_renamed.g.vcf.gz"), path("merged_renamed.g.vcf.gz.tbi")
 
	script:
	"""
	bcftools annotate \
		--rename-chrs $params.chr_map \
		--write-index="tbi" \
		--output-type z \
		--output "merged_renamed.g.vcf.gz"
	"""
}

process VCF2PHY {

	publishDir { "${params.pubdir}/tree/" }, mode: 'copy'

	input:
	tuple val(ID), path(auto_vcf), path(auto_vcf_idx)

	output:
	path "*phy", emit: phy

	script:
	"""
	vcf2phylip.py \
		--min-samples-locus 2 \
		--outgroup $params.outgroup \
		--input $auto_vcf
	# vcf2phylip can also mark an outgroup; this will be written as the first in the file, which is recognized by IQTree!
	"""
}

// docker pull staphb/iqtree2:2.2.2.7
process MLtree_SNP {

	label 'medium'

        container 'docker://staphb/iqtree2:2.2.2.7'

        input:
        path phy

        publishDir { "${params.pubdir}/tree/" }, mode: 'copy'

        output:
        path "${phy.baseName}.iqtree", emit: iqtree
        path "${phy.baseName}.treefile", emit: newick
        path "${phy.baseName}.log", emit: log

        script:
        """
        iqtree2 \
                -T AUTO \
                -s $phy \
                -m TEST \
                -bb 1000

        # default root is first taxon in alignment; this will be put py vcf2phylip if the outgroup is specified there. useful!
	# automatically select best model after AIC; this is default behaviour
        # GTR+ASC should be the correct one for SNP data; I could also define this.
        """
}

process TREE_DRAW {

	label 'small'

	container 'docker://choelter/r-base_ape:v5.8'

	input:
	path newick

        publishDir { "${params.pubdir}/tree/" }, mode: 'copy'

	output:
	path "tree_mqc.pdf"

        shell:
        '''
        #!/usr/bin/env Rscript
        library(ape)

	tree <- read.tree(text = !{newick})

	pdf(file="tree_mqc.pdf")
	plot(tree, "radial", use.edge.length=TRUE, show.tip.label=TRUE)
	dev.off() 
	'''
}

// MutliQC

process MULTIQC {

	label 'mini'

        container 'https://depot.galaxyproject.org/singularity/multiqc:1.13a--pyhdfd78af_1'

        input:
	path('*')

        publishDir { "${params.pubdir}/multiqc/" }, mode: 'copy'
        
        output:
	path "multiqc_report.html", emit: multiqc_report
	path "multiqc_data", emit: multiqc_data

	script:
        """
        multiqc .
        """
}

// Workflow call

workflow {
	// create empty channels to fill later
	multiqc_files   = Channel.from([])

	// these should implicitly be value channels ...
	fa_ch = Channel.fromPath(params.refg_fa, checkIfExists: true)
	gtf_ch = Channel.fromPath(params.refg_gtf, checkIfExists: true)
	cdna_ch = Channel.fromPath(params.refg_gtf, checkIfExists: true)

        Channel
                .fromPath(params.samplesheet, checkIfExists:true)
                .splitCsv(sep:",", header: true)
                .map{
                        row -> tuple(row.SAMPLE_ID, row.RUN_ID, row.R1, row.R2)
                }
                .set{ch_reads}

	FAIDX(fa_ch)	
	DICT(fa_ch)

	// making value ch
	ref_fa_ch = FAIDX.out.ref_fa.collect()
	ref_fai_ch = FAIDX.out.ref_fai.collect()
	ref_dict_ch = DICT.out.ref_dict.collect()

	Mapping(ref_fa_ch, ch_reads)
	BamProc(ref_fa_ch, ref_fai_ch, ref_dict_ch, Mapping.out.clean_fq, Mapping.out.bam_aln)

	bam_ch = BamProc.out.RMDUP_BAM
	bam_ch.view()
	// take all calculated COV values, and concat the min cov parameter - then find the minimum:
	//min_cov_ch = bam_ch.map({ it[1] })
	//	.concat(Channel.value(params.min_cov))
	//	.min()
	//	.view()

	//bam_ch = bam_ch.filter { it[1].toFloat() > min_cov_ch.val }

	// troubleshooting
	//println "Min Cov: ${params.min_cov}"

	//bam_ch_filtered = bam_ch.filter { it[1].toFloat() > params.min_cov }
	//bam_ch_filtered.view()
	//bam_ch_filtered.ifEmpty { warn "ERROR: There were no BAM files that met your criteria for coverage: ${params.min_cov}" }

	// I should make a report sheet and note it therre if there were no files left!

	// orient myself on minimal coverage; and if the lowest mean coverage in the dataset is above that, take that as minimum to be downsampled to.
	// do a bool here to decide if we want to do downsampling
	//if ( params.downsampling ){
	//	DOWNSAMPLING(bam_ch)
	//}

	// remove the COV value from channel, as it is no longer needed.
	bam_ch_reduced = bam_ch.map { [it[0], it[2], it[3]] }

	SNPcall(ref_fa_ch, ref_fai_ch, ref_dict_ch, bam_ch_reduced)
	// SNP_anno here
	//SNPanno(ref_fa_ch, ref_fai_ch, gtf_ch, cdna_ch, SNPcall.out.SNP)

	// tree
	VCF_MERGE(SNPcall.out.SNP.collect())

	// PCA
	PCA(VCF_MERGE.out)

	VCF2PHY(VCF_MERGE.out)
	MLtree_SNP(SNPcall.out.SNP)
	TREE_DRAW(MLtree_SNP.out.newick)

	// write csv file with IDs and vcf file paths - to read in with second part of pipeline - i.e. admixx
	SNPcall.out.SNP.collectFile() { id, vcf, idx -> ["${params.pubdir}/vcf_files.csv", "${id},${vcf},${idx}\n"] }

	// QC
	// QC here, for input files
	// rather than giving all channels individually, mix everything together with .mix() step by step.
	// this way, empty channels do not hold up the execution of multiqc! - hopefully, needs to be tested!
	//QualiWorker(
	//	multiqc_files.mix(
	//		ch_reads
	//		Mapping.out.clean_fq,
	//		Mapping.out.bam_aln,
	//		bam_ch
	//	)
	//)

	// add to MULTIQC channel
	multiqc_files.mix(
		Mapping.out.fastp_qc,
		BamProc.out.RMDUP_PICARD_METRICS,
                SNPcall.out.QC,
                PCA.out,
                TREE_DRAW.out

	)
	MULTIQC(multiqc_files)
}

