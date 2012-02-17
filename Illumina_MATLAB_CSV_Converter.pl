#! /usr/bin/perl -w

#Convert from CSV given to us by Illumina to CSV file for MATLAB
#Make a couple files - make a file for the targetid names - with two things in it, the illumina Tar#getID and my GeneID
#Then make a file with the sample names(column identifiers)
#Finally, make a file with all the values of interest

#Define Subroutines:
sub tumnorm {
  $rep=$sample_break[1];
  if ($rep =~ m/T/) {
    $rep=~s/T//;
    print SAMPLE $j*2,",$rep,-1\n";
  }
  if ($rep =~ m/N/) {
    $rep=~s/N//;
    print SAMPLE ($j*2)+1, ",$rep,-1\n";
  }
}




#Get Input file - make output filename
$input_file = shift;
$target_file = $input_file;
$target_file =~ s/.csv/_mat_targer.csv/;
$sample_file = $input_file;
$sample_file =~ s/.csv/_mat_sample.csv/;
$data_file = $input_file;
$data_file =~ s/.csv/_mat_data.csv/;


#Create Convert Table
$naming_file = shift;

open(CONVERT, "< $naming_file");
while (<CONVERT>) {
  #Split input
  @convert_fields = split /,/;
  @convert_min = ($convert_fields[0],$convert_fields[10]);
  push @convert_table, [@convert_min];
}

close(CONVERT);



open(CSVFILE, "< $input_file");
open(TARGET, ">", $target_file);
open(SAMPLE, ">", $sample_file);
open(DATA, ">", $data_file);

#Take input from first line - header
$header=<CSVFILE>;
@header_field = split(/,/,$header);
$header_size = @header_field;

#Take only the avg values, which are every fifth column starting at 15 - extract the name

for ($i=15; $i<$header_size; $i=$i+5) {
  $sample_name = $header_field[$i];
  #Remove trailing garbage
  $sample_name =~ s/@\S+//;
  print SAMPLE "$sample_name,";
  #Also set a class:
  #1 - Controls
  #2 - Breast Tumor
  #3 - Breast Normal
  #4 - Colon Tumor
  #5 - Colon Normal
  #6 - Lung Tumor
  #7 - Lung Normal
  #8 - Ovary Tumor
  #9 - Ovary Normal
  #10 - Wilms Tumor
  #11 - Wilms Normal
  #12 - Cross Control
  #And set the sample number and % methylation

  @sample_break = split(/_/,$sample_name);
  if ($sample_name =~ m/CIDR/) {
    #Check for replicate
    if ($sample_name =~ m/0_/) {
      print SAMPLE "1,$sample_break[2],$sample_break[1]\n";
    } else {
      print SAMPLE "1,0,$sample_break[1]\n";
    }
    
  }
  if ($sample_name =~ m/Half/) {
    $rep=$sample_break[1];
    $rep=~s/HalfM//;
    print SAMPLE "1,$rep,50\n";
  }
  if ($sample_name =~ m/Cross/) {
    print SAMPLE "12,$sample_break[2],$sample_break[3]\n";
  }
  if ($sample_name =~ m/Breast/) {
    $j=1;
    &tumnorm;
  }
  if ($sample_name =~ m/Colon/) {
    $j=2;
    &tumnorm;
  }
  if ($sample_name =~ m/Lung/) {
    $j=3;
    &tumnorm;
  }
  if ($sample_name =~ m/Ovary/) {
    $j=4;
    &tumnorm;
  }
  if ($sample_name =~ m/Wilms/) {
    $j=5;
    &tumnorm;
  }
  
}



close(SAMPLE);

#
while (<CSVFILE>) {
  s/\r\n//;
  @fields = split /,/;

  
  #Description - to do this, I need to match the ID to the match file
  foreach $aref (@convert_table) {
    if (@$aref[0] =~ m/$fields[0]/) {
	print TARGET "@$aref[1],$fields[0]\n";
	last;
      } 
    }
  #Target_ID

  #Data
  for ($i=15; $i<$header_size; $i=$i+5) {
    #AVG
    print DATA "$fields[$i],";
    #Color 1
    print DATA "$fields[$i+1],";
    #Color 2
    print DATA "$fields[$i+2],";
    #Num Beads
    print DATA "$fields[$i+3],";
    #Bead STD
    print DATA "$fields[$i+4],";

  }

  print DATA "\n";

}


close DATA;
close CSVFILE;
close TARGET;









	
    



