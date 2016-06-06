#!/usr/bin/env bash
#$ -j y
#$ -pe smp 8
#$ -V

source ${HOME}/.bashrc

module load novoalign/3.01.02
module load bioinformatics/samtools/0.1.19
module load picard-tools/2.2.2
module load bedtools/2.16.2
module load ucsc_tools/13_05_2016
module load general/python/2.7.10

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

    cp "${ref}" "${TMPDIR}/${ref##*/}"
    REFERENCE="${TMPDIR}/${ref##*/}"
    echo "New reference location: ${REFERENCE}"

else
    echo "Could not continue as no reference was provided"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit
fi

ANNOTATION=""
echo "Annotation Source = ${ann}"
if [ "${ann}" != "" ]; then

    cp "${ann}" "${TMPDIR}/${ann##*/}"
    ANNOTATION="${TMPDIR}/${ann##*/}"
    echo "New annotation location ${ANNOTATION}"

else
    echo "Could not continue as no annotation set was provided."

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

mkdir -p "${SAMPLE_OUTPUT_PATH}"

############################### ALIGNMENT TO THE REFERENCE GENOME #######################################
printf "Started Alignment on $date\n\n"
novoalign -F STDFQ -f "${READ_1}" -d "${REFERENCE}" -o SAM -o SoftClip --Q2Off -k -a -g 65 -x 7 -c 8 2> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_alignment_metrics.txt" | samtools view -bS - > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
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

############################# BAM TO BED FORMAT #########################################################
printf "Started converting to BED format on `date`\n\n"
bamToBed -i "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed"
printf "Finished converting to BED format on `date`\n\n"
#########################################################################################################

total_final_reads=`cat "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" | wc -l`
printf "\n\n\ntotal_final_reads = $total_final_reads\n\n\n" >> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_alignment_metrics.txt"

################################################ RUNNING MACS2 to make BigWig files #####################
printf "Started running MACS2 on `date`\n\n"
macs2 callpeak -t "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed" -f BED -g hs --keep-dup=all --outdir "${SAMPLE_OUTPUT_PATH}/macs2" -n "${SAMPLE_NAME}" -B --SPMR --nomodel --extsize=200 -q 0.05
printf "Finished running MACS2 on `date`\n\n"

printf "Started making genome browser supported files on `date`\n\n"
bedClip "${SAMPLE_OUTPUT_PATH}/macs2/${SAMPLE_NAME}_treat_pileup.bdg" "${ANNOTATION}" "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_treat_temp.bdg"
bedGraphToBigWig "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_treat_temp.bdg" "${ANNOTATION}" "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bw"
printf "Finished making genome browser supported files on `date`\n\n"
##########################################################################################################

printf "Tidying up on `date`\n\n"

rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam"
rm "${SAMPLE_OUTPUT_PATH}/macs2/${SAMPLE_NAME}_treat_pileup.bdg"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_treat_temp.bdg"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bed"
rm -rf "${SAMPLE_OUTPUT_PATH}/macs2"

# ================================================== CORE SAMPLE LOOP =================================================

rm "${TMPDIR}/${ref##*/}"
rm "${TMPDIR}/${ann##*/}"