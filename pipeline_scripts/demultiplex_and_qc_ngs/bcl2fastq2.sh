#!/usr/bin/env bash
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#S -V

# TODO: Parse CSV
SUBMISSION_DIRECTORY=""
OUTPUT_DIRECTORY=""

# Output directory == working directory as the file structure will need fixing in module 2
EXECUTION_VARIABLES=""

if [[ ${create_fastq_index} ]]; then
    EXECUTION_VARIABLES+="--create-fastq-for-index-reads \ "
fi

if [[ ${ignore_missing_bcl} ]]; then
    EXECUTION_VARIABLES+="--ignore-missing-bcls \ "
fi

if [[ ${ignore_missing_filter} ]]; then
    EXECUTION_VARIABLES+="--ignore-missing-filter \ "
fi

if [[ ${ignore_missing_positions} ]]; then
    EXECUTION_VARIABLES+="--ignore-missing-positions \ "
fi

if [[ ${with_failed_reads} ]]; then
    EXECUTION_VARIABLES+="--with-failed-reads \ "
fi

if [[ ${write-rev-comp} ]]; then
    EXECUTION_VARIABLES+="--write-fastq-reverse-complement \ "
fi

if [[ ${no_compression} ]]; then
    EXECUTION_VARIABLES+="--no-bgzf-compression \ "
fi

if [[ ${no_lane_splitting} ]]; then
    EXECUTION_VARIABLES+="--no-lane-splitting \ "
fi

if [[ ${find_adapters_using_sliding_window} ]]; then
    EXECUTION_VARIABLES+="--find-adapters-with-sliding-window \ "
fi

if [[ ${tiles}} ]]; then
    EXECUTION_VARIABLES+="--tiles ${tiles} \ "
fi

if [[ ${base_mask}} ]]; then
    # Split semicolons and create multiple

    EXECUTION_VARIABLES+="--base-mask ${base_mask} \ "
fi

EXECUTION_VARIABLES+="
    --runfolder-dir ${WORKING_DIRECTORY} \
    --output_dir ${WORKING_DIRECTORY} \
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

module load bcl2fastq2/2.17.1.14
bcl2fastq2 "${EXECUTION_VARIABLES}"
