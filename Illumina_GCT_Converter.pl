#! /usr/bin/perl -w

#Convert from CSV given to us by Illumina to GCT file for GenePattern
#Also make CLS File


#Get Input file - make output filename
$input_file = shift;
$gct_file = $input_file;
$gct_file =~ s/.csv/.gct/;
$cls_file = $input_file;
$cls_file =~ s/.csv/.cls/;

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
open(OUTPUT, ">", $gct_file);
open(CLSFILE, ">", $cls_file);

#Take input from first line - header
$header=<CSVFILE>;
@header_field = split(/,/,$header);
$header_size = @header_field;
push @rower, "NAME";
push @rower, "Description";
#Take only the avg values, which are every fifth column starting at 15
for ($i=15; $i<$header_size; $i=$i+5) {
  push @rower, $header_field[$i];
}
$num_col = @rower-2;

push @fully, [@rower];

#Setup CLS file - first #samples, then #classes, then 1
print CLSFILE "$num_col 3 1\n";
#Put in class names
print CLSFILE "# tumor normal other \n";
#Remove the Name and Description lines
$garb=shift(@rower);
$garb=shift(@rower);
#Check for Tumor or Normal - assign classes
foreach (@rower) {
  if (m/_T/) {
    print CLSFILE "0 ";
  } else {
    if (m/_N/) {
      print CLSFILE "1 ";
    } else {
      print CLSFILE "2 ";
    }
  }
  
}
print CLSFILE "\n";
close CLSFILE;
	     


#
while (<CSVFILE>) {
  $#rower = -1;
  @fields = split /,/;

  
  #Description - to do this, I need to match the ID to the match file
  foreach $aref (@convert_table) {
    if (@$aref[0] =~ m/$fields[0]/) {
	push @rower, @$aref[1];
	last;
      } 
    }
  #Target_ID
  push @rower, $fields[0];  
  #Data
  for ($i=15; $i<$header_size; $i=$i+5) {
    push @rower, $fields[$i];
  }

  push @fully, [@rower];
}

close CSVFILE;





$num_row = @fully-1;

print OUTPUT "#1.2\n";
print OUTPUT "$num_row\t$num_col\n";
for $aref ( @fully ) {
  foreach (@$aref) {
    print OUTPUT "$_\t";
  }
  print OUTPUT "\n";
}

close OUTPUT









	
    



