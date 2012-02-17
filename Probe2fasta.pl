#! /usr/bin/perl -w

#Load in the old probe file, and convert it to fasta format.

#Get Input data file - make output filename
$input_file = shift;
$target_file = $input_file;
$target_file =~ s/.csv/.fa/;

#Load Gene info

open(CONVERT, "<", "$input_file");

open(FASTA, ">", "$target_file");

$/ = "\r";

$garb=<CONVERT>;
$garb=0;

while (<CONVERT>) {
  #Split input
  @convert_fields = split /,/;
  print FASTA ">$convert_fields[0]\n";
  print FASTA "$convert_fields[1]\n";
}

close(FASTA);
close(CONVERT);



