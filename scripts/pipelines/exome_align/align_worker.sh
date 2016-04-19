#!/usr/bin/env bash
#$ -j y
#$ -pe threaded 8

module load bioinformatics/novalign
module load bioinformatics/samtools/0.1.19
module load bioinformatics/picardtools
module load bioinformatics/bedtools

OLDIFS="${IFS}"
date=`date`

# Debug
echo "Array ID: ${SGE_TASK_ID}"
echo "Data File: ${DATA_FILE}"

# =========================================== BUILD OUR EXECUTION VARIABLES! ==========================================
echo "Beginning runtime arguments parsing..."

REFERENCE=""
echo "reference_genome = ${ref}"
if [ "${ref}" != "" ]; then

    cp "${ref}" "${TMPDIR}"
    REFERENCE="${TMPDIR}/${ref##*/}"

else
    echo "Could not continue as no reference was provided"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit
fi
# =================================== DONE BUILDING OUR EXECUTION VARIABLES! ==========================================

# ================================================== CORE SAMPLE LOOP =================================================
echo "Beginning alignment module"

# Get the line we are interested in
IFS=','
while read SAMPLE_NAME SAMPLE_INPUT_PATH SAMPLE_OUTPUT_PATH EXTRA
do

    READ_1=""
    READ_2=""

    # Expectation of nested pipeline datasets within sample dir (current layout of biocompute - may change, see the two errors below)
    for d in "${SAMPLE_INPUT_PATH}/*/"; do

        # We are looking for a specific file type
        for f in d/*_R1.fastq; do
            if [ "${READ_1}" ]; then
                echo "More that one R1 was identified. Exome alignment cannot determine which you wish to align. This is a programming error and indicative of a current flaw in Biocompute that will be addressed asap"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else

                READ_1="${f}"
            fi
        done

        # We are looking for a specific file type
        for f in d/*_R2.fastq; do
            if [ "${READ_2}" ]; then
                echo "More that one R2 was identified. Exome alignment cannot determine which you wish to align. This is a programming error and indicative of a current flaw in Biocompute that will be addressed asap"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                READ_2="${f}"
            fi
        done
    done

    if [ -z "${READ_1}" ]; then
        echo "Could not identifiy the primary read fastq"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

        exit
    else
        echo "Read 1 was identified as: ${READ_1}"
    fi

    if [ "${READ_2}" ]; then
        echo "Read 2 was identified as: ${READ_2}"
    fi

    ############################### ALIGNMENT TO THE REFERENCE GENOME #######################################
    printf "Started Alignment on $date\n\n"
    novoalign -F STDFQ -f "${READ_1}" "${READ_2}" -d "${REFERENCE}" -o SAM -o SoftClip --Q2Off -k -a -g 65 -x 7 -c 8 2> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_alignment_metrics" | samtools view -bS - > "${SAMPLE_INPUT_PATH}/${SAMPLE_NAME}.bam"
    printf "Finished Alignment on $date\n\n"
    #########################################################################################################

    ############################# SORTING THE BAM FILE  #####################################################
    printf "Started sorting the BAM file on $date\n\n"
    java -Xmx30g -jar picard.jar SortSam INPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam" OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
    printf "Finished sorting the BAM file on $date\n\n"
    #########################################################################################################

    ############################# REMOVE PCR DUPLICATES  ####################################################
    printf "Started removing PCR duplicates on $date\n\n"
    java -Xmx30g -jar picard.jar MarkDuplicates I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam" M="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_duplicate_metrics" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
    printf "Finished removing PCR duplicates on $date\n\n"
    #########################################################################################################

    ############################# INDEX BAM FILES  ##########################################################
    printf "Started indexing BAM file on $date\n\n"
    samtools index "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam"
    printf "Finished indexing BAM file on $date\n\n"
    #########################################################################################################

    ############################# COVERAGE METRICS  V2 ######################################################

    # coverage calculations
    samtools view -bq 20 -F 1796 "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" | bamToBed -i stdin > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed"
    coverageBed -hist -a "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -b bed_files/targets.bed > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_targets_cov.bed" &
    coverageBed -a "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -b bed_files/baits.bed > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits_cov.bed" &
    coverageBed -a "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -b bed_files/baits_plus_150bp.bed > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits150_cov.bed" &
    wait

    #efficiency of capture
    reads=`awk 'END {OFS = "\t"; print NR}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" `
    mapped_to_target_reads=`awk '{SUM += $4} END {OFS = "\t";print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits_cov.bed"`
    percent1=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads' / '$reads') * 100)}'`
    mapped_to_target_reads_plus_150=`awk '{SUM += $4} END {OFS = "\t";print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits150_cov.bed"`
    percent2=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads_plus_150' / '$reads') * 100)}'`

    #coverage
    grep all "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_targets_cov.bed" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"
    meancov=`awk '{if ($2>=1) (SUM += $2*$5)} END {printf ("%0.2f", SUM)}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`

    #completeness of coverage
    cov1xpc=`awk '{if ($2>=1) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov5xpc=`awk '{if ($2>=5) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov10xpc=`awk '{if ($2>=10) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov20xpc=`awk '{if ($2>=20) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov0x=`awk '{if ($2>=0) (SUM += $3)} END {print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov1x=`awk '{if ($2>=1) (SUM += $3)} END {print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov5x=`awk '{if ($2>=5) (SUM += $3)} END {print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov10x=`awk '{if ($2>=10) (SUM += $3)} END {print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`
    cov20x=`awk '{if ($2>=20) (SUM += $3)} END {print SUM}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"`

    #report generation
    printf "total_reads\t"$reads"\n" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "mapped_to_target_reads\t"$mapped_to_target_reads"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "percentage\t"$percent1"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "mapped_to_target_reads_plus_150bp\t"$mapped_to_target_reads_plus_150"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "percentage\t"$percent2"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "mean_coverage\t"$meancov"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "accessible_target_bases\t"$cov0x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "accessible_target_bases_1x\t"$cov1x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "percentage_1x\t"$cov1xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "accessible_target_bases_5x\t"$cov5x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "percentage_5x\t"$cov5xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "accessible_target_bases_10x\t"$cov10x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "percentage_10x\t"$cov10xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "target_bases_20x\t"$cov20x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    printf "percentage_20x\t"$cov20xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics"
    #####################################################################################################

    printf "Tidying up on $date\n\n"

    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam"
    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed"
    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_targets_cov.bed"
    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits_cov.bed"
    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits150_cov.bed"
    rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"

done < `awk "NR==${SGE_TASK_ID}" "${DATA_FILE}"`
# ================================================== CORE SAMPLE LOOP =================================================

rm "${TMPDR}/${ref##*/}"