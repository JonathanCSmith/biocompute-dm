#!/usr/bin/env bash
#$ -j y
#$ -pe smp 8
#$ -V

source ${HOME}/.bashrc

module load novoalign/3.01.02
module load bioinformatics/samtools/0.1.19
module load picard-tools/2.0.1
module load bedtools/2.16.2
module load gatk/3.4-0

OLDIFS="${IFS}"
date=`date`

# Debug
echo "Array ID: ${SGE_TASK_ID}"
echo "Data File: ${DATA_FILE}"

# =========================================== BUILD OUR EXECUTION VARIABLES! ==========================================
echo "Beginning runtime arguments parsing..."

REF=""
echo "reference_genome = ${ref}"
if [ "${ref}" != "" ]; then

    cp -a "${ref}" "${TMPDIR}/${ref##*/}"
    REF="${TMPDIR}/${ref##*/}"
    echo "New reference location: ${REF}"

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
IFS=',' read SAMPLE_NAME SAMPLE_INPUT_PATH SAMPLE_OUTPUT_PATH EXTRA < <(sed -n ${SGE_TASK_ID}p ${DATA_FILE})

READ_1=""
READ_2=""
REGEX=".*\\*.*"

# We are looking for a specific file type
for f in ${SAMPLE_INPUT_PATH}/*_R1_001.fastq; do
    if [[ "${f}" =~ $REGEX ]]; then
        echo "Ignoring: ${f} as it is likely a false positive"

    elif [ "${READ_1}" ]; then
        echo "More that one R1 was identified. Exome alignment cannot determine which you wish to align."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

        exit

    else
        echo "Identified ${f}"
        READ_1="${f}"
    fi
done

# We are looking for a specific file type
for f in ${SAMPLE_INPUT_PATH}/*_R1_001.fastq.gz; do

    if [[ "${f}" =~ $REGEX ]]; then
        echo "Ignoring: ${f} as it is likely a false positive"

    elif [ "${READ_1}" ]; then
        echo "More that one R1 was identified. Exome alignment cannot determine which you wish to align"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

        exit

    else
        echo "Identified ${f}"
        READ_1="${f}"
    fi
done

# We are looking for a specific file type
for f in ${SAMPLE_INPUT_PATH}/*_R2_001.fastq; do
    if [[ "${f}" =~ $REGEX ]]; then
        echo "Ignoring: ${f} as it is likely a false positive"

    elif [ "${READ_2}" ]; then
        echo "More that one R2 was identified. Exome alignment cannot determine which you wish to align"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

        exit

    else
        echo "Identified ${f}"
        READ_2="${f}"
    fi
done

# We are looking for a specific file type
for f in ${SAMPLE_INPUT_PATH}/*_R2_001.fastq.gz; do
    if [[ "${f}" =~ $REGEX ]]; then
        echo "Ignoring: ${f} as it is likely a false positive"

    elif [ "${READ_2}" ]; then
        echo "More that one R2 was identified. Exome alignment cannot determine which you wish to align"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

        exit

    else
        echo "Identified ${f}"
        READ_2="${f}"
    fi
done

if [ -z "${READ_1}" ]; then
    echo "Could not identifiy the primary read fastq"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit

elif [[ "${READ_1}" =~ $REGEX ]]; then
    exit

else
    echo "Read 1 was identified as: ${READ_1}"
fi

if [ "${READ_2}" ]; then
    echo "Read 2 was identified as: ${READ_2}"

    if [[ "${READ_2}" =~ $REGEX ]]; then
        echo "Read 2 was likely malformed as it contained a regex pointer. Read 2 will be discarded"
        READ_2 = ""
    fi
fi

mkdir -p "${SAMPLE_OUTPUT_PATH}"


############################### ALIGNMENT TO THE REFERENCE GENOME #######################################
printf "Started Alignment on $date\n\n"
novoalign -F STDFQ -f "${READ_1}" "${READ_2}" -d $REF/novoindex/novoindex -o SAM -o SoftClip --Q2Off -k -a -g 65 -x 7 -c 8 2> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_alignment_metrics.txt" | samtools view -bS - > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
printf "Finished Alignment on $date\n\n"
#########################################################################################################

############################# SORTING THE BAM FILE  #####################################################
printf "Started sorting the BAM file on $date\n\n"
java -Xmx30g -jar ${PICARD}/picard.jar SortSam INPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam" OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
printf "Finished sorting the BAM file on $date\n\n"
#########################################################################################################

############################# REMOVE PCR DUPLICATES  ####################################################
printf "Started removing PCR duplicates on $date\n\n"
java -Xmx30g -jar ${PICARD}/picard.jar MarkDuplicates I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" M="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_duplicate_metrics.txt" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
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
coverageBed -hist -a "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -b $REF/bedfiles/wes_targets.bed > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_targets_cov.bed" &
coverageBed -a "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -b $REF/bedfiles/wes_baits.bed > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits_cov.bed" &
coverageBed -a "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -b $REF/bedfiles/wes_baits_padded.bed > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits150_cov.bed" &
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
printf "total_reads\t"$reads"\n" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "mapped_to_target_reads\t"$mapped_to_target_reads"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "percentage\t"$percent1"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "mapped_to_target_reads_plus_150bp\t"$mapped_to_target_reads_plus_150"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "percentage\t"$percent2"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "mean_coverage\t"$meancov"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "accessible_target_bases\t"$cov0x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "accessible_target_bases_1x\t"$cov1x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "percentage_1x\t"$cov1xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "accessible_target_bases_5x\t"$cov5x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "percentage_5x\t"$cov5xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "accessible_target_bases_10x\t"$cov10x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "percentage_10x\t"$cov10xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "target_bases_20x\t"$cov20x"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
printf "percentage_20x\t"$cov20xpc"\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage_metrics.txt"
#####################################################################################################


############################# Post Alignment QC Metrics ##########################################################

java -Xmx30g -jar ${PICARD}/picard.jar AddOrReplaceReadGroups I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final_mod.bam" RGLB="Library1" RGPL="illumina" RGPU="HISEQ3000" RGSM=$sample

samtools index "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final_mod.bam"


java -Xmx30g -jar ${GATK}/GenomeAnalysisTK.jar -T DiagnoseTargets -R $REF/refgenome/refgenome.fa -I "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final_mod.bam" -L $REF/bedfiles/wes_targets.bed -missing "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_missing_intervals.txt" -o "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_DiagnoseTargets.txt"


java -Xmx30g -jar ${GATK}/GenomeAnalysisTK.jar -T DepthOfCoverage -R $REF/refgenome/refgenome.fa -o "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_DepthOfCoverage.txt" -I "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final_mod.bam" -L $REF/bedfiles/wes_targets.bed


##################################################################################################################

printf "Tidying up on $date\n\n"

rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_targets_cov.bed"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits_cov.bed"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_baits150_cov.bed"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_coverage.hist"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final_mod.bam"
# ================================================== CORE SAMPLE LOOP =================================================

rm "${TMPDIR}/${ref##*/}"
