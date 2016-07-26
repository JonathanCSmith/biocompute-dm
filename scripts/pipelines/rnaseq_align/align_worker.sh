#!/usr/bin/env bash
#$ -j y
#$ -pe smp 8
#$ -V

source ${HOME}/.bashrc

module load novoalign/3.01.02
module load bioinformatics/samtools/0.1.19
module load picard-tools/2.2.2
module load bedtools/2.16.2

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

    cp -a "${ref}" "${TMPDIR}/${ref##*/}"
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
        echo "More that one R1 was identified. RNAseq alignment cannot determine which you wish to align."

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
        echo "More that one R1 was identified. RNAseq alignment cannot determine which you wish to align"

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
        echo "More that one R2 was identified. RNAseq alignment cannot determine which you wish to align"

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
        echo "More that one R2 was identified. RNAseq alignment cannot determine which you wish to align"

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
printf "Started Alignment on `date`  \n\n"
hisat2 -q --trim5 12 --known-splicesite-infile $REF/hisat2index/splice_sites --dta-cufflinks -p 8 --seed 3 -x $REF/hisat2index/hisat2index -1 "${READ_1}" -2 "${READ_2}" -S "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.sam" 2> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_alignment_metrics.txt"
printf "Finished Alignment on `date`  \n\n"
#########################################################################################################

############################# CONVERTING SAM TO BAM FORMAT ##############################################
printf "Started converting to BAM format on `date`  \n\n"
samtools view -bS "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.sam" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
printf "Finished converting to BAM format on `date`  \n\n"
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


############################# OTHER QUALITY METRICS  ##########################################################

printf "Started colelcting metrics on `date`  \n\n"

java -Xmx30G -jar ${PICARD}/picard.jar CollectRnaSeqMetrics REF_FLAT=$REF/refflat/transcript_all STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND CHART_OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_chart.pdf" I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_RnaSeqMetrics.txt"

java -Xmx30G -jar ${PICARD}/picard.jar CollectRnaSeqMetrics REF_FLAT=$REF/refflat/transcript_upto4Kb STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND CHART_OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_chart_upto4Kb.pdf" I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_RnaSeqMetrics_upto4Kb.txt"

java -Xmx30G -jar ${PICARD}/picard.jar CollectRnaSeqMetrics REF_FLAT=$REF/refflat/transcript_4To8Kb STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND CHART_OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_chart_4To8Kb.pdf" I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_RnaSeqMetrics_4To8Kb.txt"

java -Xmx30G -jar ${PICARD}/picard.jar CollectRnaSeqMetrics REF_FLAT=$REF/refflat/transcript_8Kb STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND CHART_OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_chart_8Kb.pdf" I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_RnaSeqMetrics_8Kb.txt"

printf "Finished collecting metrics on `date`  \n\n"

###############################################################################################################


printf "Tidying up on $date\n\n"

rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.sam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam"

# ================================================== CORE SAMPLE LOOP =================================================

rm "${TMPDIR}/${ref##*/}"
