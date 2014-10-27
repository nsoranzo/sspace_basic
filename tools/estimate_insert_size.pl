###################################################################################################################
#Marten Boetzer BaseClear B.v. 14-07-2011                                                                         #
#SSPACE perl subscript samToTab_multi.pl                                                                          #
#This script;                                                                                                     #
#  -Estimates median insert size by mapping paired-reads on contigs                                               #
#  It goes through each contig and maps both reads, if a pair is mapped,                                          #
#  the orientation and insert size is estimated.                                                                  #
#  If sufficient pairs (given by the user) are found, the median insert size is                                   #
#  estimated, as well as a file with the distribution is generated which can be                                   #
#  used to visualize the insert size distribution.                                                                #
#                                                                                                                 #
#  To run this script;                                                                                            #
#  perl estimate_insert_size.pl <contigfile> <readfile1> <readfile2> <number_of_pairs> <orientation_of_pairs>     #
                                                                                                                  #
#  Output is the median insert size and a file with distribution of the insert size. Also, number of pairs for    #
#  each found orientation (FR, RF, FF and RR) are given.                                                          #
###################################################################################################################

use File::Path;
use strict;
my $contigfile = $ARGV[0];
my $fileA = $ARGV[1];
my $fileB = $ARGV[2];
my $numpairs = $ARGV[3];
my $orientation = $ARGV[4];

die "ERROR: Can't find contig file: $contigfile -- fatal\n" if(! -e $contigfile);
die "ERROR: Can't find read file: $fileA -- fatal\n" if(! -e $fileA);
die "ERROR: Can't find read file: $fileB -- fatal\n" if(! -e $fileB);
if($numpairs eq ''){
  print "WARNING: No number of pairs are given, using 10000 pairs instead\n";
  $numpairs = 10000;
}
if($orientation eq ''){
  print "WARNING: No orientation of the pairs is given, using orientation FR instead\n";
  $orientation = "FR";
}
die "ERROR: You've inserted $numpairs, which does not seem to be an valid number. Exiting.\n" if(!($numpairs>0) || !($numpairs =~ /^\d+$/));
die "ERROR: Orientation must have length of 2 characters and should contain one of the following; FR, FF, FR or RF. You've inserted orientation of $orientation ...Exiting.\n" if(!(length($orientation) == 2) || !($orientation =~ /[FR][FR]/));

print "\n";
my $paircount = 0;
my ($direction, $insertsize);
mkpath('bowtieoutput');
open (CONT, $contigfile) || die "Can't open contig file $contigfile\n";

my ($seq,$name, $maxctg, $maxseq, $maxname)=("","",0,"","");
my $contignum = 0;
CONTIG:
while (<CONT>) {
  chomp;
  $seq.=$_ if(eof(CONT));
  if (/\>(\S+)/ || eof(CONT)){
    if($seq ne ""){         
       $contignum++;
       if(length($seq) > $maxctg){
         $maxctg = length($seq);
         $maxseq = $seq;
         $maxname = $name;
       }
       if(eof(CONT)){
         $seq = $maxseq;
         $name = $maxname;
       }
       if(eof(CONT)){
         print "now at contig $name = size".length($seq)."\n";
         open (BOWCONT, ">bowtieoutput/bowtie_input.fa");
         print BOWCONT ">$name\n$seq\n";
         close BOWCONT;
         ($paircount) = &mapWithBowtie($contignum,"bowtieoutput/bowtie_input.fa", $fileA, $fileB);
         last CONTIG if($paircount>=$numpairs);
       }

       $name = "";
       $seq = "";
    }
    $name = $1;
  }
  else {
     $seq .= $_;
  }
}

foreach my $d (keys %$direction){
  print "direction $d is found $direction->{$d} times\n";
}
my ($median_ins,$record) = (0,0);
my $median_bin = int($paircount/2);
open (CSV, ">distribution.txt") || die "Can't open distribution.txt for writing -- fatal";
foreach my $is (sort {$a<=>$b} keys %$insertsize){
  for(my $i=0;$i<$insertsize->{$is};$i++){
    $record++;
    $median_ins = $is if($record >= $median_bin && $median_ins == 0);
  }
  print CSV "$is\t$insertsize->{$is}\n";
}

print "\nmedian = $median_ins\n\nSee the distribution in file 'distribution.txt'\n";


sub mapWithBowtie{
  my ($fname,$contig, $fileA, $fileB) = @_;
  my $bowtieout = "contig$fname.bowtieIndex";
  system("bowtie-build $contig bowtieoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?"; # returns exit status values

  my $fastq = 0;
  open(TEST, "< $fileA");
  $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);

  open(FILEA, "< $fileA");
  open(FILEB, "< $fileB");

  my $count=0;
  open (BOWIN, ">bowtieoutput/bowtiein.$fname.fa") || die "Can't write to single file bowtieoutput/bowtiein.$fname.fa-- fatal\n";
  while(<FILEA>) {
    <FILEB>;
    $count++;
    my $seq1 = <FILEA>;
    chomp $seq1;
    my $seq2 = <FILEB>;
    chomp $seq2;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);
    
    print BOWIN ">read$count\n$seq1>read$count\n$seq2";
    if($count > $numpairs){
      close BOWIN;
      open(IN, "bowtie -p 1 -v 0 -m 1 bowtieoutput/$bowtieout --suppress 6,7 -f bowtieoutput/bowtiein.$fname.fa --quiet|") || die "Can't open bowtie output -- fatal\n";
      my ($prevread, $prevline);
      while(my $line = <IN>){
        my @t1 = split(/\t/,$line);
        if($prevread eq $t1[0]){
          $paircount++;
          my @t2 = split(/\t/,$prevline);
          my ($start1, $start2, $end1,$end2);

          if($t1[1] eq "+"){
            $end1 = $t1[3] + length($t1[4]);
            $start1 = $t1[3];
          }
          else{
            $start1 = $t1[3] + length($t1[4]);
            $end1 = $t1[3];
          }
          if($t2[1] eq "+"){
            $end2 = $t2[3] + length($t2[4]);
            $start2 = $t2[3];
          }
          else{
            $start2 = $t2[3] + length($t2[4]);
            $end2 = $t2[3];
          }
          my ($dir1, $dir2);
          $dir1 = "F" if($start1 < $end1);
          $dir1 = "R" if($start1 > $end1);
          $dir2 = "F" if($start2 < $end2);
          $dir2 = "R" if($start2 > $end2);
          $direction->{"$dir1$dir2"}++ if($start1 < $start2);
          $direction->{"$dir2$dir1"}++ if($start2 < $start1);
          my $diff = abs($start2-$start1);
          if($orientation eq "$dir1$dir2" || $orientation eq "$dir2$dir1"){
            $insertsize->{$diff}++;
          }
          return $paircount if($paircount >= $numpairs);
        }
        $prevread = $t1[0];
        $prevline = $line;
     }

      close BOWIN;
      open (BOWIN, "bowtieoutput/bowtiein.$fname.fa") || die "Can't write to single file bowtieoutput/bowtiein.$name.fa-- fatal\n";
    }
  }
  print "count = $paircount\n";
  return $paircount;
}

###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}
