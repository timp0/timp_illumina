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
	$ori_namey = $1;
	if ($ori_namey =~ /(\S+)_/) {
	    $namey=$1;
	    #Quotemeta puts backslahes before control variables
	    $ori_namey= quotemeta($ori_namey);
	    $linein =~ s/$ori_namey/$namey/;
	} else {
	    $namey = $ori_namey;
	}
	open(OUTPUT, "> $namey.gb");
	
    }
    
    if ($linein !~ m/(ORGANISM|SOURCE|KEYWORDS|VERSION|\s\.\s\n)/) {
	print OUTPUT $linein;
    }
    if (/\/\//) {
	close OUTPUT;
    }
}


