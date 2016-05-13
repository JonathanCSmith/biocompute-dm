#!/usr/bin/env bash
#$ -j y
#$ -pe smp 8
#$ -V

source ${HOME}/.bashrc

module load hisat2/2.0.3
module load bioinformatics/samtools/0.1.19
module load picard-tools/2.2.2
module load bedtools/2.16.2
module load ucsc_tools/13_05_2016

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

SPLICE=""
echo "splice sites = ${splice}"
if [ "${splice}" != "" ]; then

    cp "${splice}" "${TMPDIR}/${splice##*/}"
    SPLICE="${TMPDIR}/${splice##*/}"
    echo "New reference location: ${SPLICE}"

else
    echo "Could not continue as no splice information was provided"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit
fi

ANNOTATIONS=""
echo "splice sites = ${ann}"
if [ "${ann}" != "" ]; then

    cp "${ann}" "${TMPDIR}/${ann##*/}"
    ANNOTATIONS="${TMPDIR}/${ann##*/}"
    echo "New reference location: ${ANNOTATIONS}"

else
    echo "Could not continue as no annotations were provided"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit
fi

GENOME=""
echo "splice sites = ${gen}"
if [ "${gen}" != "" ]; then

    cp "${gen}" "${TMPDIR}/${gen##*/}"
    GENOME="${TMPDIR}/${gen##*/}"
    echo "New reference location: ${GENOME}"

else
    echo "Could not continue as no genome was provided"

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

# We are looking for a specific file type
for f in ${d}/*_R1.fastq; do
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

if [ -z "${READ_1}" ]; then
    echo "Could not identifiy the primary read fastq"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit

elif [[ "${READ_1}" =~ ".**.*" ]]; then
    exit

else
    echo "Read 1 was identified as: ${READ_1}"
fi


############################### ALIGNMENT TO THE REFERENCE GENOME #######################################
printf "Started Alignment on `date`  \n\n"
hisat2 -q --trim5 12 --known-splicesite-infile "${SPLICE}" --novel-splicesite-outfile "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_novel_splice_sites" --dta-cufflinks -p 8 --seed 3 -x "${REFERENCE}" -1 $2 -2 $3 -S "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.sam" 2> "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_alignment_metrics.txt"
printf "Finished Alignment on `date`  \n\n"
#########################################################################################################


############################# CONVERTING SAM TO BAM FORMAT ##############################################
printf "Started converting to BAM format on `date`  \n\n"
samtools view -bS "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.sam" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
printf "Finished converting to BAM format on `date`  \n\n"
#########################################################################################################


############################# SORTING THE BAM FILE  #####################################################
printf "Started sorting the BAM file on `date`  \n\n"
java -Xmx30G -jar ${PICARD}/picard.jar SortSam INPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam" OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
printf "Finished sorting the BAM file on `date`  \n\n"
#########################################################################################################


############################# REMOVE PCR DUPLICATES  ####################################################
printf "Started removing PCR duplicates on `date`  \n\n"
java -Xmx30G -jar ${PICARD}/picard.jar MarkDuplicates I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" M="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_duplicate_metrics.txt" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
printf "Finished removing PCR duplicates on `date`  \n\n"
#########################################################################################################


############################# INDEX BAM FILES  ##########################################################
printf "Started indexing BAM file on `date`  \n\n"
samtools index "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam"
printf "Finished indexing BAM file on `date`  \n\n"
#########################################################################################################


############################# OTHER QUALITY METRICS  ##########################################################
printf "Started colelcting metrics on `date`  \n\n"
java -Xmx30G -jar ${PICARD}/picard.jar CollectRnaSeqMetrics REF_FLAT="${ANNOTATIONS}" STRAND_SPECIFICITY=NONE CHART_OUTPUT="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_chart.pdf" I="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" O="${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_RnaSeqMetrics.txt"
printf "Finished collecting metrics on `date`  \n\n"
###############################################################################################################


############################# GENERATE BIGWIG FILES  ####################################################
printf "Started making BIGWIG file on `date`  \n\n"
genomeCoverageBed -split -bg -ibam "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" -g "${GENOME}" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bdg"

##### Remove the ERCC lines from the bedGraph file #####
grep -v ERCC "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bdg" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_mod.bdg"


##### Normalize the bedGraph to library size #####
total_size=`samtools view "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bam" | wc -l`
awk '{FS=OFS="\t"; print $1,$2,$3,($4/'$total_size')*1000000}' "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_mod.bdg" > "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bdg"


##### Converting bedGraph to bigWig #####
bedGraphToBigWig "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bdg" "${GENOME}" "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bw"
printf "Finished making BIGWIG file on `date`  \n\n"
#########################################################################################################

printf "Tidying up on `date`  \n\n"

rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.sam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_sorted.bam"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}.bdg"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_mod.bdg"
rm "${SAMPLE_OUTPUT_PATH}/${SAMPLE_NAME}_final.bdg"

# ================================================== CORE SAMPLE LOOP =================================================

rm "${TMPDIR}/${ref##*/}"
rm "${TMPDIR}/${splice##*/}"
rm "${TMPDIR}/${ann##*/}"
rm "${TMPDIR}/${gen##*/}"