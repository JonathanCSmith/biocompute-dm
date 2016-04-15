#!/usr/bin/env bash
#$ -j y
#$ -pe threaded 8

module load bioinformatics/novalign
module load bioinformatics/samtools/0.1.19
module load bioinformatics/picardtools
module load bioinformatics/bedtools

OLDIFS="${IFS}"
date=`date`

# =========================================== BUILD OUR EXECUTION VARIABLES! ==========================================
echo "Beginning runtime arguments parsing..."

REFERENCE=""
echo "reference_genome = ${ref}"
if [ "${ref}" != "" ]; then
    REFERENCE="${ref}"
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

# Loop through our sample information and process
IFS=','
while read SAMPLE_NAME SAMPLE_INPUT_PATH SAMPLE_OUTPUT_PATH EXTRA
do

    # Path to the sample directory
    sample="${SAMPLE_INPUT_PATH}"

done < ${SAMPLE_CSV}
# ================================================== CORE SAMPLE LOOP =================================================

sample="sample_"$1

mkdir $sample

############################### ALIGNMENT TO THE REFERENCE GENOME #######################################

printf "Started Alignment on $date\n\n"

$NOVOALIGN/novoalign -F STDFQ -f $2 $3 -d $NOVOINDEX/TBW_ref_genome.index -o SAM -o SoftClip --Q2Off -k -a -g 65 -x 7 -c 8 2> $sample/$sample"_alignment_metrics" | $SAMTOOLS/samtools view -bS - > $sample/sample.bam

printf "Finished Alignment on $date\n\n"

#########################################################################################################


############################# SORTING THE BAM FILE  #####################################################

printf "Started sorting the BAM file on $date\n\n"

java -Xmx30g -jar $PICARDTOOLS/picard.jar SortSam INPUT=$sample/sample.bam OUTPUT=$sample/sample_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

printf "Finished sorting the BAM file on $date\n\n"

#########################################################################################################


############################# REMOVE PCR DUPLICATES  ####################################################

printf "Started removing PCR duplicates on $date\n\n"

java -Xmx30g -jar $PICARDTOOLS/picard.jar MarkDuplicates I=$sample/sample_sorted.bam O=$sample/$sample"_final.bam" M=$sample/$sample"_duplicate_metrics" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

printf "Finished removing PCR duplicates on $date\n\n"

#########################################################################################################


############################# INDEX BAM FILES  ##########################################################

printf "Started indexing BAM file on $date\n\n"

$SAMTOOLS/samtools index $sample/$sample"_final.bam"

printf "Finished indexing BAM file on $date\n\n"

#########################################################################################################


############################# COVERAGE METRICS  V1 ##########################################################

#printf "Started generating coverage metrics on $date\n\n"

#java -Xmx30g -jar $PICARDTOOLS/picard.jar CalculateHsMetrics BI=bed_files/baits_test.bed TI=bed_files/targets_test.bed INPUT=$sample/$sample"_final.bam" OUTPUT=$sample/$sample"_coverage_metrics"

#printf "Finished generating coverage metrics on $date\n\n"

#########################################################################################################



############################# COVERAGE METRICS  V2 ##########################################################

# coverage calculations
$SAMTOOLS/samtools view -bq 20 -F 1796 $sample/$sample"_final.bam" | $BEDTOOLS/bamToBed -i stdin > $sample/$sample"_final.bed"
$BEDTOOLS/coverageBed -hist -a $sample/$sample"_final.bed" -b bed_files/targets.bed > $sample/$sample"_targets_cov.bed" &
$BEDTOOLS/coverageBed -a $sample/$sample"_final.bed" -b bed_files/baits.bed > $sample/$sample"_baits_cov.bed" &
$BEDTOOLS/coverageBed -a $sample/$sample"_final.bed" -b bed_files/baits_plus_150bp.bed > $sample/$sample"_baits150_cov.bed" &
wait

#efficiency of capture
reads=`awk 'END {OFS = "\t"; print NR}' $sample/$sample"_final.bed" `
mapped_to_target_reads=`awk '{SUM += $4} END {OFS = "\t";print SUM}' $sample/$sample"_baits_cov.bed"`
percent1=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads' / '$reads') * 100)}'`
mapped_to_target_reads_plus_150=`awk '{SUM += $4} END {OFS = "\t";print SUM}' $sample/$sample"_baits150_cov.bed"`
percent2=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads_plus_150' / '$reads') * 100)}'`

#coverage
grep all $sample/$sample"_targets_cov.bed" > $sample/$sample"_coverage.hist"
meancov=`awk '{if ($2>=1) (SUM += $2*$5)} END {printf ("%0.2f", SUM)}' $sample/$sample"_coverage.hist"`

#completeness of coverage
cov1xpc=`awk '{if ($2>=1) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' $sample/$sample"_coverage.hist"`
cov5xpc=`awk '{if ($2>=5) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' $sample/$sample"_coverage.hist"`
cov10xpc=`awk '{if ($2>=10) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' $sample/$sample"_coverage.hist"`
cov20xpc=`awk '{if ($2>=20) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' $sample/$sample"_coverage.hist"`
cov0x=`awk '{if ($2>=0) (SUM += $3)} END {print SUM}' $sample/$sample"_coverage.hist"`
cov1x=`awk '{if ($2>=1) (SUM += $3)} END {print SUM}' $sample/$sample"_coverage.hist"`
cov5x=`awk '{if ($2>=5) (SUM += $3)} END {print SUM}' $sample/$sample"_coverage.hist"`
cov10x=`awk '{if ($2>=10) (SUM += $3)} END {print SUM}' $sample/$sample"_coverage.hist"`
cov20x=`awk '{if ($2>=20) (SUM += $3)} END {print SUM}' $sample/$sample"_coverage.hist"`

#report generation
printf "total_reads\t"$reads"\n" > $sample/$sample"_coverage_metrics"
printf "mapped_to_target_reads\t"$mapped_to_target_reads"\n" >> $sample/$sample"_coverage_metrics"
printf "percentage\t"$percent1"\n" >> $sample/$sample"_coverage_metrics"
printf "mapped_to_target_reads_plus_150bp\t"$mapped_to_target_reads_plus_150"\n" >> $sample/$sample"_coverage_metrics"
printf "percentage\t"$percent2"\n" >> $sample/$sample"_coverage_metrics"
printf "mean_coverage\t"$meancov"\n" >> $sample/$sample"_coverage_metrics"
printf "accessible_target_bases\t"$cov0x"\n" >> $sample/$sample"_coverage_metrics"
printf "accessible_target_bases_1x\t"$cov1x"\n" >> $sample/$sample"_coverage_metrics"
printf "percentage_1x\t"$cov1xpc"\n" >> $sample/$sample"_coverage_metrics"
printf "accessible_target_bases_5x\t"$cov5x"\n" >> $sample/$sample"_coverage_metrics"
printf "percentage_5x\t"$cov5xpc"\n" >> $sample/$sample"_coverage_metrics"
printf "accessible_target_bases_10x\t"$cov10x"\n" >> $sample/$sample"_coverage_metrics"
printf "percentage_10x\t"$cov10xpc"\n" >> $sample/$sample"_coverage_metrics"
printf "target_bases_20x\t"$cov20x"\n" >> $sample/$sample"_coverage_metrics"
printf "percentage_20x\t"$cov20xpc"\n" >> $sample/$sample"_coverage_metrics"


#########################################################################################################


printf "Tidying up on $date\n\n"

rm $sample/sample.bam
rm $sample/sample_sorted.bam
rm $sample/$sample"_final.bed"
rm $sample/$sample"_targets_cov.bed"
rm $sample/$sample"_baits_cov.bed"
rm $sample/$sample"_baits150_cov.bed"
rm $sample/$sample"_coverage.hist"

