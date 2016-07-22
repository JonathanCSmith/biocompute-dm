#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

module load general/R/3.2.1

echo "Executing R script for plot generation"
Rscript "${PIPELINE_SOURCE}"/plot_ChIP_qc_metrics.R "./merged_qc_table" "./pipeline_output/"
echo "Plotting complete"

mv "./merged_qc_table" "./pipeline_output/merged_qc_table.txt"