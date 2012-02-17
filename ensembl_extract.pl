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


$q=1;
$k=0;


#Take input from first line
$garb=<>;
#print $garb,"\n";

#Take input
while (<>) {
    #Split input
    @fields = split /,/;
    #/Now read in all this stuff
    $id_name = $fields[5];
    $chromey = $fields[0];
    $begin = $fields[1]-$span;
    $finish = $fields[2]+$span;
    $delta_m = $fields[3];
    $fdr = $fields[4];
    
    $relation = $fields[6];
    $TSS_dist = $fields[7];
    $CGI = $fields[8];
    $CGI_dist = $fields[9];
    #$index = $11;

    for ($i=0; $i<=(($#fields-10)/2); $i++) {
	$probe[$i]=$fields[2*$i+11];
	$deltam[$i]=$fields[2*$i+12];
    }
    print "$id_name $probe[0] \n";

    
    
	
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
			    -accession_number => sprintf("%03d",$q),
			    -desc => "$CGI:$CGI_dist Chromosome $chromey $begin-$finish DeltaM: $delta_m FDR: $fdr $relation to $id_name, $TSS_dist bp");
    
    
    #Define annotation for CHARM Region - this was the original start and end given in the line
    $feat = new Bio::SeqFeature::Generic(-start => $span,
					 -end => ($sender->length)-$span,
					 -primary_tag => 'CHARM_Region',
					 -tag => {note => 'Geneious name: CHARM Region'});
    $sender->add_SeqFeature($feat);
    
    
#Probes are 50bp long - Define Annotation for Probe(s)
    
    for ($i=0; $i<=(($#probes); $i++) {
	$feat = new Bio::SeqFeature::Generic(-start => ($probe[$i]-$begin),
					     -end => ($probe[$i]-$begin+50),
					     -primary_tag => "CHARM_Probe_$i",
					     -tag => {DeltaM => $deltam[$i], note => 'Geneious name: CHARM Probe #$i'});
	$sender->add_SeqFeature($feat);
	 }
		
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
	    
	    $il_start=$curprobe;
	    $il_end=$curprobe+50;
	    
	    print "$CGI,$chromey,$il_start,$il_end,36.1,Homo sapiens,$id_name.$i\n";

	    $sequency = $sender->subseq(($curprobe-$begin),($curprobe-$begin+50));
	    
	    while ($sequency =~ m/CG/g) {
		#Off by 3 - one for zero correction - two for the two bp
		#Add one on either side - get 4 bp with CG in the middle
		#$loco = pos($sequency)+$curprobe-3-1;
		#$second = $loco+1+2;
		#print "$chromey,$loco,$second,$id_name.$j\n";
		$j++;
	    }
	$k=$k+$j;
        }
	
	
	
	
	$io = Bio::SeqIO->new(-format => "genbank", file => ">$id_name.gb");
	$io->write_seq($sender);
	
	#print "$id_name has $cpgs CpGs.\n";
	
	$q++;
	

    }
}

print "For a total of $k CpGs.\n";











	
    



