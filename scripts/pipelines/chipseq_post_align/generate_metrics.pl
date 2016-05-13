open OUT, "> merged_qc_table";

print OUT "sampleNames,totalReads,alignedReads,uniquelyAlignReads,uniquelyAlignReadsMinusdup,percentDuplicates\n";

@samples = glob("sample_*");

foreach $sample (@samples)
{
	chomp $sample;

	$sample_name = $1 if $sample =~ /sample_(.+)/;

	print OUT "$sample_name,";

	$aln_file = $sample . "_alignment_metrics.txt";
	$dup_file = $sample . "_duplicate_metrics.txt";

	open ALN, "$sample/$aln_file";

	while ($line = <ALN>)
	{
		chomp $line;

		$tot = $1 if $line =~ /Read Sequences:    (\d+)/;
		$aln = $1 if $line =~ /   Aligned:    (\d+)/;
		$ualn = $1 if $line =~ /Unique Alignment:    (\d+)/;
		$final = $1 if $line =~ /total_final_reads = (\d+)/;
	}

	open DUP, "$sample/$dup_file";

	while ($line1 = <DUP>)
	{
		chomp $line1;
		$metric = $line1 if $line1 =~ /^Unknown Library/;
		@info = split("\t", $metric);
		$dup = $info[7] * 100;
	}

	print OUT "$tot,$aln,$ualn,$final,$dup\n";
}
