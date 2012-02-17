#! /usr/bin/perl -w

# First need to install EnsEMBL packages which are found at 
# www.ensembl.org/info/docs/api/
# May also need to install some DBI and DBD/msql libraries - 
# just google for the problem 
#

use Bio::Perl;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use GD;
use GD::Graph;
use GD::Graph::area;
use GD::Text;
use Bio::EnsEMBL::Registry;
use warnings;

#How much to span on either side of CHARM region. . . 
$span=500;



#Take input from first argument
while (<>) {
    if (/chr(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+),(\S+)/)
    {
	$id_name = $6;
	$chromey = $1;
	$begin = $2-$span;
	$finish = $3+$span;
	$delta_m = $4;
	$fdr = $5;
	$consist = $7;
	$largest = $8;
	
	
	
	
#Connect to ENSEMBL
	my $registry = 'Bio::EnsEMBL::Registry';
	
	$registry->load_registry_from_db(
	    -host => 'ensembldb.ensembl.org',
	    -user => 'anonymous'
	    );
	
	
	
	my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
	
	
	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chromey, $begin, $finish);
	
	
	
	
	
	
	$sender = Bio::Seq->new(-seq => $slice->seq,
				-display_id => $id_name,
				-desc => "Chromosome $chromey $begin-$finish DeltaM: $delta_m FDR: $fdr");
	
	$feat = new Bio::SeqFeature::Generic(-start => $span,
					     -end => ($sender->length)-$span,
					     -primary_tag => 'CHARM_Region',
					     -tag => {note => 'Geneious name: CHARM Region'});

	$sender->add_SeqFeature($feat);

#Probes are 50bp long
	$feat = new Bio::SeqFeature::Generic(-start => ($consist-$begin),
					     -end => ($consist-$begin+50),
					     -primary_tag => 'CHARM_TProbe',
					     -tag => {note => 'Geneious name: CHARM Consistant Probe'});

	$sender->add_SeqFeature($feat);

	$feat = new Bio::SeqFeature::Generic(-start => ($largest-$begin),
					     -end => ($largest-$begin+50),
					     -primary_tag => 'CHARM_DProbe',
					     -tag => {note => 'Geneious name: CHARM Largest Probe'});

	
	$sender->add_SeqFeature($feat);
	
	$sequency = $sender->subseq(($consist-$begin),($consist-$begin+50));
	
	#print $sequency, "\n";

	my $cpgs = ($sequency =~ s/CG/CG/g);


	$genes = $slice->get_all_Genes();
	
	
	
	
	while ( $gene = shift @{$genes} ) {
	    $feat = new Bio::SeqFeature::Generic(-start => $gene->start,
						 -end => $gene->end,
						 -strand => $gene->strand,
						 -primary_tag => 'gene',
						 -tag => {gene => $gene->external_name});
	    
	    $sender->add_SeqFeature($feat);
	}
	
	
	$io = Bio::SeqIO->new(-format => "genbank", file => ">$id_name.gb");
	$io->write_seq($sender);
	
	print "$id_name has $cpgs CpGs.\n";

    }
}











	
    



