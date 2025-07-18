#!/usr/bin/env Rscript

# This script generates a karyotype plot showing the frequency (number of samples) each genomic interval appears in.
# The input must be a BED file containing 'chrom', 'start', 'end', and 'num' columns.

# Usage: Rscript karyoplot.R <path/to/intervals.bed>
# Note: This script runs only on HPC/Linux. On other systems, it may throw: Error in options(scipen = old.scipen): invalid 'scipen'

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)        # for data analysis, manipulation, and visualization (loads dplyr, tidyr, ggplot2, readr, stringr, etc.)
  library(karyoploteR)      # to create karyotype plots for genomes
  library(GenomicRanges)    # for manipulating genomic intervals
  library(viridis)  # for color mapping
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript karyotype.R <path/to/intervals.bed>")
file_path <- args[1]
if (!file.exists(file_path)) stop("Error: Common intersections BED file not found!")

# Load BED file into a data frame, ensuring necessary columns exist
bed_data <- read.table(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!all(c("chrom", "start", "end", "num") %in% colnames(bed_data))) {
  stop("Error: BED file must contain columns: chrom, start, end, num")
}

# Convert BED data into a GRanges object
regions <- GRanges(seqnames = bed_data$chrom,
                   ranges = IRanges(start = bed_data$start, end = bed_data$end),
                   num = bed_data$num)  # Store read counts as metadata

# Save plot as a high-resolution image
jpeg(file.path(dirname(file_path), "/../../../results/host_dna_analysis/karyotype_frequency.jpg"), width=3000, height=3000, res=300)

# Modify plot parameters for better visualization
custom_params <- getDefaultPlotParams(plot.type = 1)
custom_params$data1outmargin <- custom_params$data1outmargin * 3  # Increase spacing between chromosomes
custom_params$leftmargin <- custom_params$leftmargin / 2  # Adjust left margin

# Generate karyotype plot
kp <- plotKaryotype(genome="hg38", plot.type=1, plot.params=custom_params)

# Add title and background for readability
kpAddMainTitle(kp, main="Metagenomic Read Alignment to the Human Genome Across Samples", cex=1)
kpDataBackground(kp, data.panel=1)

# Add chromosome names and base number scale
kpAddChromosomeNames(kp)

kpAddBaseNumbers(kp, tick.dist=20000000, minor.tick.dist=5000000, 
                 tick.len=2, minor.tick.len=1, cex=0.6, add.units=TRUE)

# Add y-axis ticks only if range is small
if (max(bed_data$num, na.rm = TRUE) < 8){
 # Add y-axis with manual ticks and labels
 kpAxis(kp, r1=1, ymin=0, ymax=max(bed_data$num), numticks=max(bed_data$num)+1,
        tick.len=15e5, side=2, cex=0.3, lwd=0.3)
}

# Plot coverage bars, scaling read counts appropriately
#color=viridis(max(mcols(regions)$num), option = "D", direction=-1)[mcols(regions)$num]
color = inferno(max(mcols(regions)$num))[mcols(regions)$num]
#color = magma(max(mcols(regions)$num))[mcols(regions)$num]
kpBars(kp, data=regions, y1=mcols(regions)$num / max(bed_data$num), lwd=2.5, col=color, border=color)

# Close the plotting device to save the output
dev.off()
