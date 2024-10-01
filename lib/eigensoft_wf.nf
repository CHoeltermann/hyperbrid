#!/usr/bin/env nextflow

# docker pull biocontainers/eigensoft:v7.2.1dfsg-1-deb_cv1
# docker pull mercury/eigensoft:latest

process ONE {

        container 'docker://biocontainers/eigensoft:v7.2.1dfsg-1-deb_cv1'

        input:

        publishDir { "" }, mode: 'copy'

        output:

        script:
        """
        """
}

// dataset/ vcf file should be: minor allele frequency > 5%, converted into plink format, LD filtering to remove linked sites
process SMART_PCA {

        container 'docker://biocontainers/eigensoft:v7.2.1dfsg-1-deb_cv1'

        input:
	path pedsnp

        //publishDir { "" }, mode: 'copy'

        output:


        script:
        """
	# convertf can be used to make proper input file for smartpca
	# PLNK's format PED can be used as input, for SNP data this would be .pedsnp

	smartpca \
		-snpname $pedsnp \
		-indivname \
		-genotypename \
		-numoutevec 25 \ # num of eigenvec to output
		-evecoutname "${.baseName}_eigenvectors" \
		-evaloutname "${.baseName}_eigenvalues" \
		> logfile.txt
        """
}

process PCA_R_VIZ {

        container 'docker://nanozoo/ggplot2:3.4.0--4856650'

        input:

        publishDir { "" }, mode: 'copy'

        output:

        script:
        """
	Rscript --vanilla src/pca_viz.R \
		$params.pubdir \
		"${.baseName}" \
		$evec \
		$eval
        """
}


workflow EIG {

        take:

        main:

        emit:
}

