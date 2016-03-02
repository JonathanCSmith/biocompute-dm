#!/usr/bin/env bash
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#S -V

OLDIFS=${IFS}

# ================================================= BUILD OUR IO VALUES ===============================================
echo "Beginning demultiplexing module"

# Note we are only interested in the first line as we are only expecting 1 folder!
IFS="," read NAME DATA_INPUT_DIRECTORY DATA_OUTPUT_DIRECTORY < <(sed -n 1p < "${samples}")
MODULE_OUTPUT_DIRECTORY="${module_output_directory}" # This is our output directory for display module specific stuffz

echo "Submission input directory: ${DATA_INPUT_DIRECTORY}"
echo "Submission output directory: ${DATA_OUTPUT_DIRECTORY}"
echo "Module output directory: ${MODULE_OUTPUT_DIRECTORY}"
# ========================================== FINISHED BUILDING OUR IO VALUES ==========================================


# =========================================== BUILD OUR EXECUTION VARIABLES! ==========================================
echo "Beginning runtime arguments parsing..."
EXECUTION_VARIABLES=""
if [[ ${create_fastq_index} == "False" ]]; then
    EXECUTION_VARIABLES+=" --create-fastq-for-index-reads"
fi

if [[ ${ignore_missing_bcl} == "False" ]]; then
    EXECUTION_VARIABLES+=" --ignore-missing-bcls"
fi

if [[ ${ignore_missing_filter} == "False" ]]; then
    EXECUTION_VARIABLES+=" --ignore-missing-filter"
fi

if [[ ${ignore_missing_positions} == "False" ]]; then
    EXECUTION_VARIABLES+=" --ignore-missing-positions"
fi

if [[ ${with_failed_reads} == "False" ]]; then
    EXECUTION_VARIABLES+=" --with-failed-reads"
fi

if [[ ${write-rev-comp} == "False" ]]; then
    EXECUTION_VARIABLES+=" --write-fastq-reverse-complement"
fi

if [[ ${no_compression} == "False" ]]; then
    EXECUTION_VARIABLES+=" --no-bgzf-compression"
fi

if [[ ${no_lane_splitting} == "False" ]]; then
    EXECUTION_VARIABLES+=" --no-lane-splitting"
fi

if [[ ${find_adapters_using_sliding_window} == "False" ]]; then
    EXECUTION_VARIABLES+=" --find-adapters-with-sliding-window"
fi

if [[ ${tiles}} == "False" ]]; then
    # Split semicolons and create multiple
    IFS=";" read -r -a array <<< "${tiles}"
    for field in "${array[@]}";
    do
        EXECUTION_VARIABLES+=" --tiles ${field}"
    done
fi

if [[ ${base_mask}} == "False" ]]; then
    # Split semicolons and create multiple
    IFS=";" read -r -a array <<< "${base_mask}"
    for field in "${array[@]}";
    do
        EXECUTION_VARIABLES+=" --base-mask ${field}"
    done
fi

EXECUTION_VARIABLES+="
    --input-dir ${DATA_INPUT_DIRECTORY} \
    --output_dir ${DATA_OUTPUT_DIRECTORY} \
    --runfolder-dir ${MODULE_OUTPUT_DIRECTORY} \
    --sample-sheet ${sample_sheet} \
    --loading_threads 2 \
    --demultiplexing_threads 4 \
    --processsing_threads 8 \
    --writing_threads 2 \
    --adapter_stringency ${adapter_stringency} \
    --aggregated-tiles ${aggregated_tiles} \
    --barcode-mismatches ${barcode_mismatches} \
    --mininmum-trimmed-read-length ${minimum_read_length} \
    --mask-short-adapter-reads ${masked_adapter_read_length} \
    --fastq-compression-level ${compression_level} \
"

echo "Calculated runtime arguments: ${EXECUTION_VARIABLES}"
# =================================== DONE BUILDING OUR EXECUTION VARIABLES! ==========================================

module load bioinformatics/bcl2fastq2/2.17.1.14
bcl2fastq "${EXECUTION_VARIABLES}"

IFS="${OLDIFS}"
