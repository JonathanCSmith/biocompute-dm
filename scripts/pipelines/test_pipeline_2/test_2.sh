#!/usr/bin/env bash

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

for id in `seq 1 10`;
do
    output_dir="${DATA_OUTPUT_DIRECTORY}//${id}"
    mkdir "${output_dir}"

    outfile="${output_dir}//${id}.txt"
    echo "test_data" > "${outfile}"

    echo "placed fake data at: ${outfile}"
done

IFS="${OLDIFS}"