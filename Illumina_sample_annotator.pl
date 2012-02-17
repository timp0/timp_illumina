#! /usr/bin/perl -w

#Add annotations to sample file


#Get Input file - make output filename
$input_file = shift;
$out_file = $input_file;
$out_file =~ s/.csv/_anno.csv/;


$anno_file = "Grades02.csv";

open(ANNO, "<", $anno_file);
open(SAMPLE, "<", $input_file);
open(OUT, ">", $out_file);


@anno = <ANNO>;
$anno_header = shift(@anno);
chomp($anno_header);
print $anno_header, "\n";

@head_field=split(/,/,$anno_header);

$s_header = <SAMPLE>;
chomp($s_header);
#Print header line
print OUT $s_header;
#Num annotations
$num_anno=0;
foreach (@head_field[1..$#head_field]) {
  print OUT ",$_";
  $num_anno=$num_anno+1;
}
print OUT "\n";

while (<SAMPLE>) {
  $unmatched=1;
  s/\n//;
  s/\r//;
  $full_sample=$_;
  @fields = split(/,/,$full_sample);
  $sample_id=$fields[0];
 GILLIGAN: foreach (@anno) {
    s/\n//;
    s/\r//;
    @anno_field = split /,/;
    if ($anno_field[0] =~ m/$sample_id/) {
      print OUT "$full_sample";
      foreach (@anno_field[1..$#anno_field]) {
	print OUT ",$_";
      }
      print OUT "\n";
      $unmatched=0;
      last GILLIGAN;
    }
  }
  if ($unmatched) {
    print  OUT "$full_sample";
    for ($i = 1; $i <=$num_anno; $i++) {
      print  OUT ",-9";
    }
    print  OUT "\n";
  }

}





close(ANNO);
close(SAMPLE);
close(OUT);



