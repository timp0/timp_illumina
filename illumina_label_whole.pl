#! /usr/bin/perl -w

# First need to install EnsEMBL packages which are found at 
# www.ensembl.org/info/docs/api/
# May also need to install some DBI and DBD/msql libraries - 
# just google for the problem 
#

use Bio::Perl;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use warnings;



#Take input from first line - header - then discard
$garb=<>;
$garb=0;
#print $garb,"\n";

#Take input
while (<>) {
    #Split input
    @fields = split /,/;
    #/Now read in all this stuff
    $probe_id = $fields[0];
    $seqy = $fields[1];
    $gbname = "$fields[2].gb";

    
    $il_score = $fields[3];
    $fail_code = $fields[4];


    $seq_object = read_sequence($gbname);

    $sequence = $seq_object->seq();
   

    #g means global replacement
    #This gets rid of brackets around central CG
    $seqy =~ s/\[//gi;
    $seqy =~ s/\]//gi;


    if ($sequence =~ /$seqy/g){
	$n = pos($sequence);
	} else {
	    $n = -1;
	}

    $m = $n-length($seqy)+1;


    $feat = new Bio::SeqFeature::Generic(-start => $m,
					 -end => $n,
					 -primary_tag => "Illumina_Probe",
					 -tag => {Illumina_id => $probe_id, Illumina_Score => $il_score, Fail_Code => $fail_code});

    if ($il_score != 0) {
	$seq_object->add_SeqFeature($feat);
    }


    $io = Bio::SeqIO->new(-format => "genbank", file => ">$gbname");
    $io->write_seq($seq_object);

}













	
    



