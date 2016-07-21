open OUT, "> merged_WES_qc_metrics_table";

print OUT "sampleNames,totalReads,alignedReads,uniquelyAlignReads,uniquelyAlignReadsMinusdup,";
print OUT "percentDuplicates,cov_1x,cov_5x,cov_10x,cov_20x,capture,capture150\n";

@samples = glob("sample_*");

foreach $sample (@samples)
{
	chomp $sample;
	
	$sample_name = $1 if $sample =~ /sample_(.+)/;
	
	print OUT "$sample_name,";
	
	$aln_file = $sample . "_alignment_metrics.txt";
	$dup_file = $sample . "_duplicate_metrics.txt";
	$cov_file = $sample . "_coverage_metrics.txt";
	
	open ALN, "$sample/$aln_file";
	
	while ($line = <ALN>)
	{
		chomp $line;
		
		$tot = $1 if $line =~ /Read Sequences:    (\d+)/;
		$aln = $1 if $line =~ /   Aligned:    (\d+)/;
		$ualn = $1 if $line =~ /Unique Alignment:    (\d+)/;
	}
	
	open DUP, "$sample/$dup_file";
	
	while ($line1 = <DUP>)
	{
		chomp $line1;
		push(@info, $line1) if $line1 =~ /^Unknown/;
		$dup = $info[7] * 100;
	}
	
	open COV, "$sample/$cov_file";
	
	while ($line2 = <COV>)
	{
		chomp $line2;
		
		$final = $1 if $line2 =~ /total_reads	(.+)/;
		$x1 = $1 if $line2 =~ /percentage_1x	(.+)/;
		$x5 = $1 if $line2 =~ /percentage_5x	(.+)/;
		$x10 = $1 if $line2 =~ /percentage_10x	(.+)/;
		$x20 = $1 if $line2 =~ /percentage_20x	(.+)/;
		
		$capture = $1 if $line2 =~ /capture-percentage	(.+)/;
		$capture150 = $1 if $line2 =~ /capture150-percentage	(.+)/;

	}
	
	print OUT "$tot,$aln,$ualn,$final,$dup,$x1,$x5,$x10,$x20,$capture,$capture150\n";
}
