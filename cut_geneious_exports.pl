#! /usr/bin/perl -w

# First need to install EnsEMBL packages which are found at 
# www.ensembl.org/info/docs/api/
# May also need to install some DBI and DBD/msql libraries - 
# just google for the problem 
#

use warnings;

#How much to span on either side of CHARM region. . . 




while (<>) {
    $linein = $_;
    if (/LOCUS       (\S+)/) {
	$namey = $1;
	if ($namey =~ /(\S+)_/) {
	    $namey=$1;
	} 
	open(OUTPUT, "> $namey.gb");
    }
    
    print OUTPUT $linein;
    
    if (/\/\//) {
	close OUTPUT;
    }
}


