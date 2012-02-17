#! /usr/bin/perl -w

#Find out how far each coordinate is from a UCSC island

use List::Util qw(min);

#First load the coordinates


#Get Input file - make output filename
$input_file = shift;
$isl_file = $input_file;
$isl_file =~ s/.csv/_ucscisl.csv/;


$ucsc_cpgfile = "ucsc_cpgisl.txt";

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
print ISL ",UCSC_Relation_to_island,UCSC_Dist_to_Island,UCSC_Island_Start,UCSC_Island_End,UCSC_NumCpGs,UCSC_PerGC,UCSC_obsExp\n";

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
    @island_field = split /\t/;
    $chromy = $island_field[0];
    $chromy =~ s/chr//;
    
    if ($chromy =~ m/$fields[2]/) {
      $isl_start = $island_field[1];
      $isl_end = $island_field[2];
      #Inside Island = 0
      #Outside Island = 2
      #Overlap Island = 1

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
  print ISL "$good_isl[1],$good_isl[2],$good_isl[5],$good_isl[8],$good_isl[9]\n";

}






close(CPGFILE);
close(PROBE);
close(ISL);









	
    



