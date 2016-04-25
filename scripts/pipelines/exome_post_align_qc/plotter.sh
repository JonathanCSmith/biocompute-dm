#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

module load general/R/3.2.1

echo "Executing R script for plot generation"
Rscript "${PIPELINE_SOURCE}"/NGS_postQC_plots.R "./merged_qc_table" "./"
echo "Plotting complete"