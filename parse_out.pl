#!/usr/bin/perl;
open(FILE, "<", 'stem_cell.txt');
@file = <FILE>;
close FILE;
%hash = ();
foreach $line (@file){
	chomp $line;
	@line_array = split(/\t/, $line);
	$spotid = $line_array[0];
	$lambda = $line_array[1];
	push(@{$hash{$spotid}}, $lambda);
}

open(FILE, "<", 'lung_normal.txt');
@file = <FILE>;
close FILE;
foreach $line (@file){
	chomp $line;
	@line_array = split(/\t/, $line);
	$spotid = $line_array[0];
	$lambda = $line_array[1];
	if(exists $hash{$spotid}){
		$hash{$spotid}[1] = $lambda;
	}
}

open(FILE, "<", 'lung_cancer.txt');
@file = <FILE>;
close FILE;
foreach $line (@file){
	chomp $line;
	@line_array = split(/\t/, $line);
	$spotid = $line_array[0];
	$lambda = $line_array[1];
	if(exists $hash{$spotid}){
		$hash{$spotid}[2] = $lambda;
	}
}

foreach $key (sort keys %hash){
	print $key,"\t", $hash{$key}[0], "\t",$hash{$key}[1], "\t", $hash{$key}[2], "\n";
}