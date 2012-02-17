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


#Take input from first argument


while (<>) {
    if (/(\S+),chr(\S+),(\S+),(\S+)/)
    {
	print "1: $1 2: $2 3: $3, 4: $4 \n";
    }
}

$id_name = $1;
$chromey = $2;
$begin = $3;
$finish = $4;


#Connect to ENSEMBL
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
    );



my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );


$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chromey, 153219037, 153220825);



print $slice->seq(), "\n";




$sender = Bio::Seq->new(-seq => $slice->seq,
			-display_id => $id_name,
                        -desc => 'Chromosome 7 153219037-153220825');

$feat = new Bio::SeqFeature::Generic(-start => 500,
				     -end => ($sender->length)-500,
                                     -primary_tag => 'misc_difference',
                                     -tag => {note => 'Geneious name: CHARM Probe'});

$sender->add_SeqFeature($feat);

$genes = $slice->get_all_Genes();




while ( $gene = shift @{$genes} ) {
    print $gene->stable_id(), "\n";
    print $gene->start(), "\n";
    print $gene->end(), "\n";
    print $gene->external_name, "\n";
    print $gene->strand(), "\n";

    $feat = new Bio::SeqFeature::Generic(-start => $gene->start,
					 -end => $gene->end,
					 -strand => $gene->strand,
					 -primary_tag => 'gene',
					 -tag => {gene => $gene->external_name});

    $sender->add_SeqFeature($feat);
}


$io = Bio::SeqIO->new(-format => "genbank", file => ">$id_name.gb");
$io->write_seq($sender);


    










	
    



