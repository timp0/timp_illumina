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


#Take three arguments from command line
#Filename which contains ids extracted from Geneious - best probes
$idfilename = shift;

#Filename which is a CSV given by illumina(w/ header removed) so I can
#extract full info for resumitting
$scorefilename = shift;

#Filename we are spitting results out to
$outfilename = shift;

#Read file into array
open(IDFILE, "< $idfilename");
@ids = <IDFILE>;



open(SCOREFILE, "< $scorefilename");

#Get rid of top header line
$garb=<SCOREFILE>;
$garb=0;

#Open outputfile
open(OUTPUT, "> $outfilename");

while(<SCOREFILE>) {
    $linein=$_;
    
    foreach $tester (@ids) {
	#get rid of \n at end of each id line
	chomp $tester;
	
	if ($linein =~ /$tester/) {
	    print OUTPUT $linein;
	}
    }
}

close(SCOREFILE);
close(IDFILE);
close(OUTPUT);










	
    



