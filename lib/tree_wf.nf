#!/usr/bin/env nextflow

process TREE_PREP {

        label 'mini'

        container 'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_1'

        input:
        path snp

        output:
        path "${snp.baseName}_autosomal.vcf.gz"

        script:
        """
        ### only autosomals SNPs?!
        # if not snp.vcf.gz, use this
        bcftools view $snp -Oz -o "${snp}.gz"
        bcftools index "${snp}.gz"

        bcftools view "${snp}.gz" \
                --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21 \
                --output-type z \
                --output "${snp.baseName}_autosomal.vcf.gz"
        """
}

process DISTMX {

        // nf does not need a python container!
        //container 'docker://python:slim-bookworm'

        input:
        path autosomal_snp

        output:
        path "*.phy", emit: phy

        script:
        """
        vcf2phylip.py --input $autosomal_snp
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
        path "${phy.baseName}.iqtree"
        path "${phy.baseName}.treefile"
        path "${phy.baseName}.log"

        script:
        """
        # auto detect how many cores 2 use
        # automatically select best model after AIC; this is default behaviour
        # GTR+ASC should be the correct one for SNP data; I could also define this.
        iqtree2 \
		-mem "${task.memory.giga}G" \
                -T AUTO \
                -s $phy \
                -m MFP \
                -B 1000
        """
}

workflow Tree {

	take:
	snp

	main:
	TREE_PREP(SNPcall.out.SNP)
	DISTMX(TREE_PREP.out)
	MLtree_SNP(TREE_PREP.out)

	emit:
	trees = MLtree_SNP.out
}

