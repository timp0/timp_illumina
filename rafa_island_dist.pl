#! /usr/bin/perl -w

#Find out how far each coordinate is from a RAFA/HMM island

use List::Util qw(min);

#First load the coordinates


#Get Input file - make output filename
$input_file = shift;
$isl_file = $input_file;
$isl_file =~ s/.csv/_hmmisl.csv/;


$ucsc_cpgfile = "hmm_isl1.txt";

open(CPGFILE, "<", $ucsc_cpgfile);
open(PROBE, "<", $input_file);
open(ISL, ">", $isl_file);


@islands = <CPGFILE>;
$isl_header = shift(@islands);
print $isl_header, "\n";


$p_header = <PROBE>;
chomp($p_header);
#Print header line
  #Name, IlluminaName, Chromosome, Start, End, Relation_to_island, Dist_to_island, Island_Start, Island_end, NumCpgs, perGC, obsExp
print ISL $p_header;
print ISL ",HMM_Relation_to_island,HMM_Dist_to_Island,HMM_Island_Start,HMM_Island_End,HMM_NumCpGs,HMM_PerGC,HMM_obsExp\n";



while (<PROBE>) {
  s/\n//;
  s/\r//;
  $full_probe=$_;
  @fields = split(/,/,$full_probe);
  $probe_start = $fields[3];
  $probe_end = $fields[4];
  $probe_dist=-1;
  
 GILLIGAN: foreach (@islands) {
    s/\n//;
    s/\r//;
    @island_field = split /,/;
    $chromy = $island_field[1];
    $chromy =~ s/"//;
    $chromy =~ s/chr//;
    
    if ($chromy =~ m/$fields[2]/) {
      $isl_start = $island_field[2];
      $isl_end = $island_field[3];
      
      if (($probe_start > $isl_start) & ($probe_end < $isl_end)) {
	$probe_loc=0;
	$probe_dist=0;
	@good_isl=@island_field;
	last GILLIGAN;
      }
      elsif ((($probe_start > $isl_start) & ($probe_start < $isl_end)) | (($probe_end > $isl_start) & ($probe_end < $isl_end))) {
	$probe_loc=1;
	$probe_dist=0;
	@good_isl=@island_field;
	last GILLIGAN;
      }
      else {
	$probe_loc=2;
	$isl_dist=min(abs($probe_end-$isl_start), abs($probe_start-$isl_end));
	if (($probe_dist<0) | ($probe_dist > $isl_dist)) {
	  $probe_dist=$isl_dist;
	  @good_isl=@island_field;
	}
      }	
    }
  }
  #Output stuff
  #CSV file that looks like:
  #Name, IlluminaName, Chromosome, Start, End, Relation_to_island, Dist_to_island, Island_Start, Island_end, NumCpgs, perGC, obsExp

  print ISL "$full_probe,";
  print ISL "$probe_loc,$probe_dist,";
  print ISL "$good_isl[2],$good_isl[3],$good_isl[5],$good_isl[7],$good_isl[8]\n";

}






close(CPGFILE);
close(PROBE);
close(ISL);









	
    



