#!/usr/bin/env nextflow

// docker pull biocontainers/eigensoft:v7.2.1dfsg-1-deb_cv1
// docker pull mercury/eigensoft:latest
// docker://biocontainers/eigensoft:v7.2.1dfsg-1-deb_cv1
// docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1

// with plink 1.9
// according to the organism used, a chr number differing from human can be specified by --autosome-sum or --chr-set

// convert to plink format
// bed, bim, fam


process PLINK_CONV {

        container 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'

        input:
	tuple path(all_vcf), path(all_vcf_idx)

        output:
	tuple val(ID), path("*bed"), path("*fam"), path("*ped"), path("*bim")

        script:
        """
	# plink --vcf genotyped.vcf.gz --make-bed --allow-extra-chr --vcf-idspace-to '_' --out testsss
	plink \
		--vcf $all_vcf \
		--make-bed \
		--allow-extra-chr \
		--vcf-idspace-to '_' \
		--out "${vcf.baseName}"
        """
}

process PLINK_PCA {

        container 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'

        input:
	tuple path(bed), path(fam), path(ped), path(bim)

        output:
	tuple path(eigenvec), path(eigenval)  

        script:
        """
	plink \
		--pca \
		--header \
		--tabs \
		--bfile "${bed.baseName}"

	#--no-fid
	#--no-parents
	#--no-sex
	#--no-pheno

        # ioutput: file.eigenvec and file.eigenval
	# visualize?

	"""
}

process DRAW_PCA {

        container 'docker://nanozoo/ggplot2:3.4.0--4856650'

        input:
	tuple path(eig_vec), path(eig_val)

        publishDir { "${params.pubdir}/pca/" }, mode: 'copy'

        output:
	// multiqc needs _mqc at the end of the filename to embed the pic. Make it small, so the HTML does not get bloated.
	path "pca_plink_mqc.pdf"

        shell:
	'''
        #!/usr/bin/env Rscript
        library("ggplot2")

        eigenval <- read.table(!{eig_val}, sep="\t", header=TRUE) 
        eigenvec <- read.table(!{eig_vec}, sep="\t", header=TRUE)

        # get tabs & header?!

        eig_perc <- round((eigenval / sum(eigenval)), 2)

        # also read in species names from input csv & add them to the plot? would be nice

        plot <- ggplot(data = eigenvec) +
                geom_point(mapping = aes(x = 'PC1', y = 'PC2', color = 'IID'), show.legend = TRUE) +
                geom_hline(yintercept = 0, linetype="dotted") +
                geom_vline(xintercept = 0, linetype="dotted") +
                labs(title = "PCA of autosomal SNPs (PLINK)", 
                        x = paste0("Principal component 1 - ", eig_perc[1,1]," %"),
                        y = paste0("Principal component 2 - ", eig_perc[2,1]," %")) +
                theme_classic()
    
        ggsave("pca_plink_mqc.pdf", plot, device="pdf", dpi=150, width=5, height=5)
        '''
}

workflow PCA {

        take:
	// tuple path("${snp.baseName}_autosomal.vcf.gz"), path("${snp.baseName}_autosomal.vcf.gz.tbi")
	vcf

        main:
	PLINK_CONV(vcf)
	PLINK_PCA(PLINK_CONV.out)
	DRAW_PCA(PLINK_PCA.out)

        emit:
	DRAW_PCA.out
}

