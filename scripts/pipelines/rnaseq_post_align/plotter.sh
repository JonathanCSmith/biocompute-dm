#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

module load general/R/3.2.1

echo "Executing R script for plot generation"
Rscript "${PIPELINE_SOURCE}"/plot_RNASeq_qc_metrics.R "./Metrics_table.csv" "./pipeline_output/"
echo "Plotting complete"

mv "./Metrics_table.csv" "./pipeline_output/Metrics_table.txt"