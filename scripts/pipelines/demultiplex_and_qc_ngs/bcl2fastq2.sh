#!/usr/bin/env bash
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#S -V

OLDIFS="${IFS}"

# ================================================= BUILD OUR IO VALUES ===============================================
echo "Beginning demultiplexing module"

# Note we are only interested in the first line as we are only expecting 1 folder!
IFS="," read NAME DATA_INPUT_DIRECTORY DATA_OUTPUT_DIRECTORY EXTRA < <(sed -n 1p < "${SAMPLE_CSV}")

echo "Submission input directory: ${DATA_INPUT_DIRECTORY}"
echo "Submission output directory: ${DATA_OUTPUT_DIRECTORY}"
echo "Extra: ${EXTRA}"
echo "Module output directory: ${MODULE_OUTPUT_DIRECTORY}"
# ========================================== FINISHED BUILDING OUR IO VALUES ==========================================


# =========================================== BUILD OUR EXECUTION VARIABLES! ==========================================
echo "Beginning runtime arguments parsing..."
EXECUTION_VARIABLES=""

# Necessary information
EXECUTION_VARIABLES+=" --output-dir ${DATA_OUTPUT_DIRECTORY}"
EXECUTION_VARIABLES+=" --runfolder-dir ${DATA_INPUT_DIRECTORY}"
EXECUTION_VARIABLES+=" --reports-dir ${MODULE_OUTPUT_DIRECTORY}"
EXECUTION_VARIABLES+=" --stats-dir ${MODULE_OUTPUT_DIRECTORY}"

echo "sample_sheet = ${sample_sheet}"
if [ "${sample_sheet}" != "SampleSheet.csv" ]; then
    EXECUTION_VARIABLES+=" --sample-sheet ${sample_sheet}"
fi

echo "create_fastq_index = ${create_fastq_index}"
if [ "${create_fastq_index}" != "False" ]; then
    EXECUTION_VARIABLES+=" --create-fastq-for-index-reads"
fi

echo "ignore_missing_bcl = ${ignore_missing_bcl}"
if [ "${ignore_missing_bcl}" != "False" ]; then
    EXECUTION_VARIABLES+=" --ignore-missing-bcls"
fi

echo "ignore_missing_filter = ${ignore_missing_filter}"
if [ "${ignore_missing_filter}" != "False" ]; then
    EXECUTION_VARIABLES+=" --ignore-missing-filter"
fi

echo "ignore_missing_positions = ${ignore_missing_positions}"
if [ "${ignore_missing_positions}" != "False" ]; then
    EXECUTION_VARIABLES+=" --ignore-missing-positions"
fi

echo "with_failed_reads = ${with_failed_reads}"
if [ "${with_failed_reads}" != "False" ]; then
    EXECUTION_VARIABLES+=" --with-failed-reads"
fi

echo "write_rev_comp = ${write_rev_comp}"
if [ "${write_rev_comp}" != "False" ]; then
    EXECUTION_VARIABLES+=" --write-fastq-reverse-complement"
fi

echo "no_compression = ${no_compression}"
if [ "${no_compression}" != "False" ]; then
    EXECUTION_VARIABLES+=" --no-bgzf-compression"
fi

echo "no_lane_splitting = ${no_lane_splitting}"
if [ "${no_lane_splitting}" != "False" ]; then
    EXECUTION_VARIABLES+=" --no-lane-splitting"
fi

echo "find_adapters_using_sliding_window = ${find_adapters_using_sliding_window}"
if [ "${find_adapters_using_sliding_window}" != "False" ]; then
    EXECUTION_VARIABLES+=" --find-adapters-with-sliding-window"
fi

echo "tiles = ${tiles}"
if [ "${tiles}" != "False" ]; then
    # Split semicolons and create multiple
    IFS=";" read -r -a array <<< "${tiles}"
    for field in "${array[@]}";
    do
        EXECUTION_VARIABLES+=" --tiles ${field}"
    done
fi

echo "base_mask = ${base_mask}"
if [ "${base_mask}" != "False" ]; then
    # Split semicolons and create multiple
    IFS=";" read -r -a array <<< "${base_mask}"
    for field in "${array[@]}";
    do
        EXECUTION_VARIABLES+=" --use-bases-mask ${field}"
    done
fi

# Has safe defaults
EXECUTION_VARIABLES+=" --adapter-stringency ${adapter_stringency} --aggregated-tiles ${aggregated_tiles} --barcode-mismatches ${barcode_mismatches} --minimum-trimmed-read-length ${minimum_read_length} --mask-short-adapter-reads ${masked_adapter_read_length} --fastq-compression-level ${compression_level}"

# Don't allow the user to change these just yet
EXECUTION_VARIABLES+=" --loading-threads 4 --demultiplexing-threads 2 --processing-threads 16 --writing-threads 4"

echo "Calculated runtime arguments: ${EXECUTION_VARIABLES}"
# =================================== DONE BUILDING OUR EXECUTION VARIABLES! ==========================================

module load bioinformatics/bcl2fastq2/2.17.1.14
bcl2fastq ${EXECUTION_VARIABLES}

IFS="${OLDIFS}"
