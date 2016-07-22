
@samples = glob("sample_*");

open OUT, "> merged_RNASeq_qc_metrics_table";
print OUT "sampleNames,totalReads,alignedReads,uniquelyAlignReads,uniquelyAlignReadsMinusdup,percentDuplicates\n";

open OALL, "> transcript_metrics_all";
open OSMALL, "> transcript_metrics_small";
open OMED, "> transcript_metrics_med";
open OLONG, "> transcript_metrics_long";
open GENE, "> genebody_metrics";

print OALL "Position\t"; print OSMALL "Position\t";
print OMED "Position\t"; print OLONG "Position\t";
print GENE "Sample\tExon\tUTR\tIntron\tIntergenic\t";

for($i=0; $i<=100; $i++)
{
	print OALL "$i\t"; print OSMALL "$i\t";
	print OMED "$i\t"; print OLONG "$i\t";
}

print OALL "\n"; print OSMALL "\n";
print OMED "\n"; print OLONG "\n";
print GENE "\n";

foreach $sample(@samples)
{
	chomp $sample;
	
	$file = $sample . "_alignment_metrics.txt";
	
	open ALIGN, "../$sample/$file";
	
	while($line = <ALIGN>)
	{
		chomp $line;
		
		$tot_reads = $1 if ($line =~ /^(\d+) reads; of these:/);
		
		$aln1 = $1 if ($line =~ /(\d+) .+ aligned concordantly exactly 1 time/);
		$aln2 = $1 if ($line =~ /(\d+) .+ aligned concordantly >1 times/);
		$aln3 = $1 if ($line =~ /(\d+) .+ aligned discordantly 1 time/);

		$aln4 = $1 if ($line =~ /(\d+) .+ aligned exactly 1 time/);
		$aln5 = $1 if ($line =~ /(\d+) .+ aligned >1 times/);
	}
	
	$aln1 = $aln1*2; $aln2 = $aln2*2; $aln3 = $aln3*2;
	
	$tot_reads = $tot_reads*2;
	$reads_aln = $aln1 + $aln2 + $aln3 + $aln4 + $aln5;
	$uniq_aln = $aln1 + $aln3 + $aln4;
	
	print OUT "$sample,$tot_reads,$reads_aln,$uniq_aln,";
	
	$dup_file = $sample . "_duplicate_metrics.txt";
	
	open DUP, "../$sample/$dup_file";
	
	while($line1 = <DUP>)
	{
		chomp $line1;
		
		if($line1 =~ /^Unknown Library/)
		{
			@array = split("\t", $line1);
			$dup = $array[7];
			
		}

	}
	
	$dup_reads = $uniq_aln*$dup;
	$uniq_aln_minus_dup = int($uniq_aln - $dup_reads);
	
	print OUT "$uniq_aln_minus_dup,$dup\n";
	
	
	############################### Transcript and GeneBody metrics #################################
	
	print OALL "$sample\t"; print OSMALL "$sample\t";
	print OMED "$sample\t"; print OLONG "$sample\t";
	print GENE "$sample\t";
	
	$all_file = $sample . "_RnaSeqMetrics.txt";
	$small_file = $sample . "_RnaSeqMetrics_upto4Kb.txt";
	$med_file = $sample . "_RnaSeqMetrics_4To8Kb.txt";
	$long_file = $sample . "_RnaSeqMetrics_8Kb.txt";
	
	open ALL, "$sample/$all_file"; @all = <ALL>;
	open SMALL, "$sample/$small_file"; @small = <SMALL>;
	open MED, "$sample/$med_file"; @med = <MED>;
	open LONG, "$sample/$long_file"; @long = <LONG>;
	
	for($i=11; $i<=111; $i++)
	{
		@array = split("\t", $all[$i]);
		chomp $array[1];
		print OALL "$array[1]\t";
		
		@array = split("\t", $small[$i]);
		chomp $array[1];
		print OSMALL "$array[1]\t";
		
		@array = split("\t", $med[$i]);
		chomp $array[1];
		print OMED "$array[1]\t";
		
		@array = split("\t", $long[$i]);
		chomp $array[1];
		print OLONG "$array[1]\t";
	}
	
	for($i=7; $i<=7; $i++)
	{
		@array = split("\t", $all[$i]);
		print GENE "$array[11]\t$array[12]\t$array[13]\t$array[14]\t";
	}
	
	print OALL "\n"; print OSMALL "\n";
	print OMED "\n"; print OLONG "\n";
	print GENE "\n";
	@array=@all=@small=@med=@long=();
	
	
}
