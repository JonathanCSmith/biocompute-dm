#!/usr/bin/env bash
echo "samples: '${samples}'"
echo "T: '${T}'"
echo "B: '${B}'"
echo "L: '${L}'"

# TODO Read samples and build an output file
INPUT="${samples}"
IFS=,
[ ! -f ${INPUT} ] && { echo "${INPUT} file not found"; exit 99; }
while read id source out
do
    outfile="${out}//${id}.txt"
    echo "test_data" > outfile
    echo "placed fake data at: ${outfile}"
done < ${INPUT}
