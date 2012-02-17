#! /usr/bin/perl -w

#Convert from CSV given to us by Illumina to CSV file for MATLAB
#Make a couple files - make a file for the targetid names - with two things in it, the illumina Tar#getID and my GeneID
#Then make a file with the sample names(column identifiers)
#Finally, make a file with all the values of interest

#Define Subroutines:
sub tumnorm {
  $other_id=$sample_break[1];
  if ($other_id =~ m/T/) {
    $other_id=~s/T//;
    $class=$j*2;
    $note=1;
  }
  if ($other_id =~ m/N/) {
    $other_id=~s/N//;
    $class=($j*2)+1;
    $note=-1;
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
$naming_file = "timp_ordered1.csv";

open(CONVERT, "< $naming_file");

#Check for near probes
$k=0;
$near_what=Inf;
$near_chrom=Q;
while (<CONVERT>) {
  #Split input
  @convert_fields = split /,/;
  if (($convert_fields[3]!~$near_chrom)|(abs($convert_fields[4]-$near_what)>10000)) {
    $k++;
    #print "$convert_fields[0],$convert_fields[10],$near_what, $near_chrom, $k\n";
    $near_what=$convert_fields[4];
    $near_chrom=$convert_fields[3];
  }
  @convert_min = ($convert_fields[0],$convert_fields[10], $k);
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

#Print Sample Header File
print SAMPLE "Sample_ID,Class,Sample_ID_2,Other_Note\n";
print DATA "Probe_ID,";

#Take only the avg values, which are every fifth column starting at 15 - extract the name

for ($i=15; $i<$header_size; $i=$i+5) {
  $sample_name = $header_field[$i];
  #Remove trailing garbage
  $sample_name =~ s/@\S+//;
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
  #13 - Other
  #4 - Adenoma - Tumor *Other_ID=-2
  #2 - DCIS (Breast) Tumor *Other_ID=-2
  #3 - DCIS (Breast) Normal *Other_ID=-2
  #3,4 - ColonSequencing Samp *Other_ID=-3
  #20 - Thyroid Tumor
  #21 - Thyroid Normal
  #ADN - OtherID=1
  #FA - OtherID=2
  #FC - OtherID=3
  #FVPTC OtherID=4
  #HA OtherID=5
  #HC OtherID=6
  #PTC OtherID=7
  #22 - Pancreas Tumor
  #23 - Pancreas Normal
  #IPMN Other_ID=1
  #NET Other_ID=2
  #LCC Other_ID=3
  #And set the sample number and % methylation
  
  $class=0;
  $note=-1;
  $other_id=-1;
  
  @sample_break = split(/_/,$sample_name);
  if ($sample_name =~ m/CIDR/) {
    $class=1;
    $note=$sample_break[1];
    #Check for replicate
    if ($sample_name =~ m/0_/) {
      $other_id=$sample_break[2];
    } else {
      $other_id=0;
    }   
  }
  if ($sample_name =~ m/Half/) {
    $class=1;
    $other_id=$sample_break[1];
    $other_id=~s/HalfM//;
    $note=50;
  }
  if ($sample_name =~ m/Cross/) {
    $class=12;
    $other_id=$sample_break[2];
    $note=$sample_break[3];
  }
  if ($sample_name =~ m/Breast/) {
    $j=1;
    &tumnorm;
  }
  if ($sample_name =~ m/DCIS/) {
    $j=1;
    &tumnorm;
    $note=2*$note;
  }
  if ($sample_name =~ m/Colon/) {
    $j=2;
    &tumnorm;
    if ($sample_name =~ m/ColonSeq/) {
      $note=3*$note;
    }
  }
  if ($sample_name =~ m/Adenoma/) {
    $j=2;
    &tumnorm;
    $note=2*$note;
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
  if ($sample_name =~ m/Thy/) {
    $j=10;
    &tumnorm;

    if ($sample_name =~ m/ThyADN/) {
      $note=1*$note;
    }
    if ($sample_name =~ m/ThyFA/) {
      $note=2*$note;
    }
    if ($sample_name =~ m/ThyFC/) {
      $note=3*$note;
    }
    if ($sample_name =~ m/ThyFVPTC/) {
      $note=4*$note;
    }
    if ($sample_name =~ m/ThyHA/) {
      $note=5*$note;
    }
    if ($sample_name =~ m/ThyHC/) {
      $note=6*$note;
    }
    if ($sample_name =~ m/ThyPTC/) {
      $note=7*$note;
    }
  }
  if ($sample_name =~ m/Pan/) {
    print "$sample_name";
    $j=11;
    &tumnorm;
    if ($sample_name =~ m/IPMN/) {
      $note=1*$note;
    }
    if ($sample_name =~ m/NET/) {
      $note=2*$note;
    }
    if ($sample_name =~ m/LCC/) {
      $note=3*$note;
    }
    print ": $note\n";
  }
    

print SAMPLE "$sample_name,$class,$other_id,$note\n";
print DATA $sample_name,"_AVG,",$sample_name,"_CY3,",$sample_name,"_CY5,",$sample_name,"_NUMBEADS,",$sample_name,"_STDERR,";


}
print DATA "\n";


close(SAMPLE);

#TARGET Header file
print TARGET "Probe_ID,Illumina_ID,Chromosome,Start_loc,Finish_loc,Good_bad,Region\n";


while (<CSVFILE>) {
  s/\r\n//;
  @fields = split /,/;

  
  #Description - to do this, I need to match the ID to the match file
  foreach $aref (@convert_table) {
    if (@$aref[0] =~ m/$fields[0]/) {
      #ProbeID, IlluminaID, chromosome, start, end(same in my case), good_probe
	print TARGET "@$aref[1],$fields[0],$fields[1],$fields[11],$fields[11], 1,@$aref[2]\n";
	print DATA "@$aref[1],";
	last;
      } 
    }

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









	
    



