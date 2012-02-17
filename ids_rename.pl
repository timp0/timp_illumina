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

#Filename which is a CSV given by illumina(w/ header removed) so I can
#extract full info for resumitting
$scorefilename = shift;

#Filename we are spitting results out to
$outfilename = shift;

open(SCOREFILE, "< $scorefilename");

#Get rid of top header line
$garb=<SCOREFILE>;
$garb=0;

#Open outputfile
open(OUTPUT, "> $outfilename");

$tester = 'NOT A GENE';
$i=0;

#Fixes the weird file annothation - very specific to this fuckup
$switch=1;


while(<SCOREFILE>) {
    $linein=$_;
    @fields = split /,/;
    #/Now read in other stuff
    $cust_anno=$fields[10];	
    if ($cust_anno !~ /$tester/) {
	if ($switch) {
	    if ($cust_anno =~ /(\S+)\.(\S+)/) {
		if ($cust_anno =~ /(\S+)\.(\S+)\.(\S+)/) {
		    $tester = "$1.$2";
		} else {
		    $tester = "$1";
		}
	    } else {
		$tester = "$cust_anno";
		$switch = 0;
	    }
	} else {
	    $tester = "$cust_anno";
	}
	$i=0;	

    }

    $newy="$tester-$i";
    $linein =~ s/$cust_anno/$newy/;

    print OUTPUT $linein;
    $i++;
}

close(SCOREFILE);
close(OUTPUT);










	
    



