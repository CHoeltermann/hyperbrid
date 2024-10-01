#!/usr/bin/env Rscript
library("ggplot2")

# take cmd line arg
args = commandArgs(trailingOnly=TRUE)

out_dir = args[1]
out_name = args[2]
eval = read.table(args[3])

evec = read.table(args[4])
#colnames(evecs) <- c("sample", paste("PC", 1:ncol(evec)-1, sep=""), "pop")

# read in smartpca data
pca_data <- data.frame()

pca_plot <- ggplot(data = pca_data, aes(x = PC1, y = PC2, label = patient))
  + geom_point(aes(colour = condition))
  + ggtitle("PCA plot")
  + xlab(paste("PC1 -", pca_var_percent[1], "%", sep=" "))
  + ylab(paste("PC2 -", pca_var_percent[2], "%", sep=" "))
  + theme_bw()

ggsave(
	filename = paste0(
		out_dir,
		exp_name,
		"/eigensoft_pca_plot.png"),
	plot = pca_plot,
	device = png
)

