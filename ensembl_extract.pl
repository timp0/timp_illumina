#! /usr/bin/perl -w

# First need to install EnsEMBL packages which are found at 
# www.ensembl.org/info/docs/api/
# May also need to install some DBI and DBD/msql libraries - 
# just google for the problem 
#

use Bio::Perl;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Bio::EnsEMBL::Registry;
use warnings;

#How much to span on either side of CHARM region. . . 
$span=500;


$i=1;
$k=0;
#Take input from first line
while (<>) {
    #print $_;
    if (/chr(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+)/)
    {
	$id_name = $6;
	$chromey = $1;
	$begin = $2-$span;
	$finish = $3+$span;
	$delta_m = $4;
	$fdr = $5;

	$relation = $7;
	$TSS_dist = $8;
	$CGI = $9;
	$CGI_dist = $10;
	$index = $11;

	$probe1 = $12;
	$deltam1 = $13;
	$probe2 = $14;
	$deltam2 = $15;
	$probe3 = $16;
	$deltam3 = $17;
	$probe4 = $18;
	$deltam4 = $19;
	$probe5 = $20;
	$deltam5 = $21;
	
	
	
	
#Connect to ENSEMBL
	my $registry = 'Bio::EnsEMBL::Registry';
	
	$registry->load_registry_from_db(
	    -host => 'ensembldb.ensembl.org',
	    -user => 'anonymous'
	    );
		
	my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
	
	#Get Specific slice based on chromosome and start and end points from input - this includes +/- span
	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chromey, $begin, $finish);
	
	
	
	#create new sequence variable from the slice
	$sender = Bio::Seq->new(-seq => $slice->seq,
				-display_id => $id_name,
				-accession_number => sprintf("%03d",$i),
				-desc => "$CGI:$CGI_dist Chromosome $chromey $begin-$finish DeltaM: $delta_m FDR: $fdr $relation to $id_name, $TSS_dist bp");
	

	#Define annotation for CHARM Region - this was the original start and end given in the line
	$feat = new Bio::SeqFeature::Generic(-start => $span,
					     -end => ($sender->length)-$span,
					     -primary_tag => 'CHARM_Region',
					     -tag => {note => 'Geneious name: CHARM Region'});
	$sender->add_SeqFeature($feat);


#Probes are 50bp long - Define Annotation for Probe(s)
	$feat = new Bio::SeqFeature::Generic(-start => ($probe1-$begin),
					     -end => ($probe1-$begin+50),
					     -primary_tag => 'CHARM Probe 1',
					     -tag => {DeltaM => $deltam1, note => 'Geneious name: CHARM Probe #1'});
	$sender->add_SeqFeature($feat);

	$feat = new Bio::SeqFeature::Generic(-start => ($probe2-$begin),
					     -end => ($probe2-$begin+50),
					     -primary_tag => 'CHARM Probe ',
					     -tag => {DeltaM => $deltam2, note => 'Geneious name: CHARM Probe #2'});
	$sender->add_SeqFeature($feat);

	$feat = new Bio::SeqFeature::Generic(-start => ($probe3-$begin),
					     -end => ($probe3-$begin+50),
					     -primary_tag => 'CHARM Probe 3',
					     -tag => {DeltaM => $deltam3, note => 'Geneious name: CHARM Probe #3'});
	$sender->add_SeqFeature($feat);

	$feat = new Bio::SeqFeature::Generic(-start => ($probe4-$begin),
					     -end => ($probe4-$begin+50),
					     -primary_tag => 'CHARM Probe 4',
					     -tag => {DeltaM => $deltam4, note => 'Geneious name: CHARM Probe #4'});
	$sender->add_SeqFeature($feat);

	$feat = new Bio::SeqFeature::Generic(-start => ($probe5-$begin),
					     -end => ($probe5-$begin+50),
					     -primary_tag => 'CHARM Probe 5',
					     -tag => {DeltaM => $deltam5, note => 'Geneious name: CHARM Probe #5'});
	$sender->add_SeqFeature($feat);




	$feat = new Bio::SeqFeature::Generic(-start => ($probe2-$begin),
					     -end => ($probe2-$begin+50),
					     -primary_tag => 'CHARM Probe 2',
					     -tag => {note => 'Geneious name: CHARM Probe #2'});
	$sender->add_SeqFeature($feat);
	
		
	#Find any genes in the slice - mark them with annotations
	$genes = $slice->get_all_Genes();
	
	while ( $gene = shift @{$genes} ) {
	    $feat = new Bio::SeqFeature::Generic(-start => $gene->start,
						 -end => $gene->end,
						 -strand => $gene->strand,
						 -primary_tag => 'gene',
						 -tag => {gene => $gene->external_name});
	    
	    $sender->add_SeqFeature($feat);
	}
	
        #Extract sequence from most consistantly different probe
	$sequency = $sender->subseq(($probe1-$begin),($probe1-$begin+50));
	
	#print $sequency, "\n";
	#Previously used to count num CpGs
	#my $cpgs = ($sequency =~ s/CG/CG/g);
	
	#stole this code:
	#print "$id_name has CpGs at: ";
	$j=0;
	for ($i=1; $i<6; $i++) {
	    if ($i==1) {
	    	$curprobe = $probe1;
	    }
	    if ($i==2) {
		$curprobe = $probe2;
	    }
	    if ($i==3) {
		$curprobe = $probe3;
	    }
	    if ($i==4) {
		$curprobe = $probe4;
	    }
	    if ($i==5) {
		$curprobe = $probe5;
	    }

	    $sequency = $sender->subseq(($curprobe-$begin),($curprobe-$begin+50));
		
	    while ($sequency =~ m/CG/g) {
		#Off by 3 - one for zero correction - two for the two bp
		#Add one on either side - get 4 bp with CG in the middle
		$loco = pos($sequency)+$curprobe-3-1;
		$second = $loco+1+2;
		print "$chromey,$loco,$second,$id_name.$j\n";
		$j++;
	    }
	$k=$k+$j;
        }
	
	
	
	
	$io = Bio::SeqIO->new(-format => "genbank", file => ">$id_name.gb");
	$io->write_seq($sender);
	
	#print "$id_name has $cpgs CpGs.\n";
	
	$i++;
	

    }
}

print "For a total of $k CpGs.\n";











	
    



