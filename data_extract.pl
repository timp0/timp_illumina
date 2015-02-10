#!/usr/bin/perl -w

##
# Extract 450k datasets to give copies
# Created on 080812 by Winston Timp

use local::lib;
use PerlIO::eol;
use Archive::Tar;
use warnings;

my $datapath=shift;

##Get the thing to match
my $check=shift;
##Location of output tarball/csv file
my $outdir=shift;
##Annofile name
my $annofile="anno.csv";
##Datafile name
my $datatar="extractdat.tgz";

##Get dir listing - go through *all* the annotation files
my @maindir=<$datapath*>;

foreach (@maindir) {
    ##Get most recent csv file
    my @annofiles=sort { -M $a <=> -M $b } <$_/*csv>;
    if (defined($annofiles[0])) {
	##print "$annofiles[0]\n";
	push @csvs, $annofiles[0];
    }    
}

##Open annotation file for our friends
open(ANNO, ">", "$outdir$annofile");
print ANNO ("Experimenter,Plate.ID,Sample.Well,Slide.ID,Array.ID,Hyb.date,Image.date,Sample.ID,Sex,Age,Race,Tissue,Status,Phenotype,Individual.ID,Source,Notes,Purification,Input.Amount\n");

##Make new tar file for data
my $tar = Archive::Tar->new();


foreach $csv (@csvs) {
    ##$csv=$csvs[8];
    ##Open this CSV file - use eol package to defeat evil mac newline (\r) problems
    open(INCSV, "<",$csv);

##Skip all header lines
    while(<INCSV>){
	##print $_;
	 if (m/Experimenter/) {	     
	     last;
	 }
    }
    
    while(<INCSV>) {
	if (m/$check/i) {
	    print ANNO $_;
	    chomp;
	    @fields=split/\,/;
	    print "$datapath, $fields[1], $fields[4]\n";
	    @needthese=<$datapath$fields[1]/$fields[4]/$fields[4]_$fields[5]*.idat>;
	    $tar->add_files( @needthese );
	}
	##last;
    }
 
    close(INCSV);
}

close(ANNO);

$tar->add_files("$outdir$annofile");

$tar->write("$outdir$datatar", COMPRESS_GZIP);
