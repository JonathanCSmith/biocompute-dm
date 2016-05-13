#!/usr/bin/perl

open OUT, "> Metrics_table.csv";

print OUT "sampleNames,totalReads,alignedReads,uniquelyAlignReads,uniquelyAlignReadsMinusdup,fractionDuplicates\n";

@samples = glob("sample_*");

foreach $sample(@samples)
{
	chomp $sample;

	$sample_name = $1 if $sample =~ /sample_(.+)/;

	print OUT "$sample_name,";

	$aln_file = $sample . "_alignment_metrics.txt";
	$dup_file = $sample . "_duplicate_metrics.txt";
	$cov_file = $sample . "_coverage_metrics.txt";

	open ALN, "$sample/$aln_file";

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

	open DUP, "$sample/$dup_file";

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
}
