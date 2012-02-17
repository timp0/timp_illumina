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



#Open directory for checking
opendir(DIR, ".");

#Look through each file
while (defined($file = readdir(DIR))) {
#Test if file is a gb file
    if ($file =~ m/\.gb/) {
	#Read file
	$seq_object = read_sequence($file);
	#Get file identifier
	$id = $seq_object->display_id();
	#Cut the _modified crap that Genious puts on
	if ($id =~ m/(\S+)_/) {
	    $id = $1;
	}
	
	#print "Display ID $id \n";
	
	#Look in the definition for the genomic start, end and chromosome
	$def = $seq_object->desc();
	if ($def =~ m/(\d+)-(\d+)/) {
	    $begin = $1;
	    #$finish = $2;
	}
	if ($def =~ m/Chromosome (\S+)/) {
	    $chromey=$1;
	}

	$i=0;
	#Go through the annotations looking for misc_features
	foreach my $feat_object ( $seq_object->all_SeqFeatures) {
	    if ($feat_object->primary_tag eq "misc_feature") {

		$gprobe_start = ($feat_object->location->start);
		$gprobe_end = ($feat_object->location->end);
		$gprobe_cpg = $gprobe_start+60+$begin;
		$first_half=$seq_object->subseq($gprobe_start, $gprobe_start+59);
		$second_half=$seq_object->subseq($gprobe_start+62, $gprobe_end);
		$seq_gprobe = "$first_half\[CG\]$second_half";
		#Find misc_feature tags - get the illumina_id tag, then output it
		for my $tag ($feat_object->get_all_tags) {             
		    if ($tag =~ m/Illumina_id/) {
			@illumina_cg = $feat_object->get_tag_values($tag);
			print "$illumina_cg[0]\n";
		    }
		} 
		#print "$id-$i,$seq_gprobe,36.1,$chromey, $gprobe_cpg,NCBI,36.1,Homo sapiens, Forward,.,., 36.1, $id-$i\n";
		$i++;

	    }
	}
    }


}
