#! /usr/bin/perl -w
# Label snp annotations - use track file I downloaded and bring them allllllll in using location file 


use Bio::Perl;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use warnings;

#Take arguments from command line
#Filename which contains list of locations and files to be annotated
$infile = shift;

#Filename which is the gb to be annotated
$snpfile = "snp131_regions";


open(SNPFILE, "< $snpfile");

while(<SNPFILE>) {
  @fieldy = split /\t/;
  # Chromosome, Start, End, RS_num, Content, Type, Hetero
  @snp_small = ($fieldy[1], $fieldy[2],$fieldy[3],$fieldy[4],$fieldy[9],$fieldy[11],$fieldy[13]);
  push @snp_table, [@snp_small];
}



open(LISTY, "< $infile");

while(<LISTY>) {
  chomp;
  @fieldy = split /\t/;
  $gbname = "$fieldy[3].gb";
  print "$gbname !!!\n";
  $seq_object= read_sequence($gbname);

  foreach $aref (@snp_table) {
    if (@$aref[0] =~ m/$fieldy[0]/) {
      if ((@$aref[1] < $fieldy[2]) & (@$aref[2] > $fieldy[1])) {
	$startloc=@$aref[1]-$fieldy[1];
	$endloc=@$aref[2]-$fieldy[1] -1;
      	$feat = new Bio::SeqFeature::Generic(-start => $startloc,
				       -end => $endloc,
				       -primary_tag => "dbSNP");
  
	$feat->add_tag_value('RS', @$aref[3]);
	$feat->add_tag_value('Value', @$aref[4]);
	$feat->add_tag_value('Type', @$aref[5]);
	$feat->add_tag_value('Hetero', @$aref[6]);
	
	$seq_object->add_SeqFeature($feat);
      }
    }
  }
  $io = Bio::SeqIO->new(-format => "genbank", file => ">$gbname");
  $io->write_seq($seq_object);

}











	
    



