#!/usr/bin/perl
open(FILE, "<", "network_cancer/GSE30784_series_matrix.txt");
@FILE = <FILE>;
close FILE;
@samples = split(/\s+/, $FILE[38]);
@status  = split(/\s+\"/, $FILE[46]);
for($i=0; $i < scalar(@status); $i++){
	$k = $samples[$i];
	$v = $status[$i];
	$k =~ s/"|\n//g;
	$v =~ s/"|\n//g;
	if($v =~ m/control$/){
		print "$k.CEL.gz\n";
	}
}