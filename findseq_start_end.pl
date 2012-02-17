#! /usr/bin/perl -w

# Get start and end for regions - put into BED file

use Bio::Perl;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use warnings;

#Get Directory Listing
@files= glob "*.gb";

foreach $which (@files) {
  $seq_object = read_sequence($which);  
  $id_name = $which;
  $id_name =~ s/.gb//;
  $definition = $seq_object->desc();

  

  @fields = split / /, $definition;
  
  @posy = split /-/, $fields[3];
  
  print "chr$fields[2]\t$posy[0]\t$posy[1]\t$id_name\n";

  



}
