########################################################################
#Marten Boetzer BaseClear B.v. 26-07-2011                              #
#SSPACE perl sam_bam2Tab.pl                                            #
#This script;                                                          #
#  -converts a .sam file to a tab file containing;                     #
#      -contig of read 1                                               #
#      -start position of read 1                                       #
#      -end position of read 1                                         #
#      -contig of read 2                                               #
#      -start position of read 2                                       #
#      -end position of read 2                                         #
#                                                                      #
#  -Sam/Bam file should contain a read pair at consecutive             #
#   lines where the first line contains the first read and             #
#   second line the second read                                        #
#   In order to have such a file, sort the sam file                    #
#   before using this script with SAMTools command:                    #
#   samtools view -uS <input.sam> | samtools sort -n - <input.sorted>  #
#                                                                      #
#  -This script requires samtools to be installed                      #
#                                                                      #
#  -Bam files should end with .bam extension                           #
#                                                                      #
#INPUT:                                                                #
#   perl sam_bam2Tab.pl <samfile> <postfixread1> <postfixread2         #
#                                                                      #
#   example:                                                           #
#   perl sam_bam2Tab.pl input.sorted.bam /1 /2 out.tab                 #
#   or                                                                 #
#   perl sam_bam2Tab.pl input.sorted.sam /1 /2 out.tab                 #
#                                                                      #
#   This means that the first read is ending with /1 while the         #
#   second read ends with /2                                           #
                                                                       #
#OUTPUT:                                                               #
#   Output of this script is saved into                                #
########################################################################

my $infile = $ARGV[0];
my $postfix1 = $ARGV[1];
my $postfix2 = $ARGV[2];
my $outfile = $ARGV[3];
die "length of postfix1 ($postfix1) has not same length of postfix2 ($postfix2). Exiting...\n" if(length($postfix1) != length($postfix2));

my $bam=($infile =~ /.bam$/)? 1:0;

if($bam){
    open(SAM, "samtools view $infile |") or die "Can't open $infile for reading -- fatal\n";
}else{
    open(SAM, "$infile") || die "Can't open $infile for reading -- fatal\n";
}
open(OUT, ">$outfile") || die "Can't open $outfile for writing -- fatal\n";

my $step = 100000;
my ($ct, $diffct, $read, $prevread, $prevline, $line);
while($line = <SAM>){
  next if($line =~ /^@/);
  ($read, undef, $chrom) = split("\t", $line);
  next if($chrom eq "*");
  if($read !~ /$postfix1$/ && $read !~ /$postfix2$/){
    warn("read $read had no suffix '$postfix1' or '$postfix2', please insert a correct suffix (e.g. '/1' and '/2')\n");
  }
  $read = substr($read,0,-(length($postfix1)));
  if($prevread eq $read){
    $pair_found++;
    my ($line1, $line2) = ($prevline,$line);
    if($prevread =~ /$postfix2$/){
      $line1 = $line;
      $line2 = $prevline;
    }
    my @arr1 = split("\t", $line1);
    my @arr2 = split("\t", $line2);

    my ($tig1,$start1,$end1, $tig2,$start2,$end2) = ($arr1[2], $arr1[3], ($arr1[3]+length($arr1[9])), $arr2[2],$arr2[3],($arr2[3]+length($arr2[9])));

    if ($arr1[1] & 16) {
      $end1 = $start1;
      $start1 = $start1 + length($arr1[9]);
    }
    if ($arr2[1] & 16) {
      $end2 = $start2;
      $start2 = $start2 + length($arr2[9]);
    }
    print OUT "$tig1\t$start1\t$end1\t$tig2\t$start2\t$end2\n";
  }
  $prevread = $read;
  $prevline = $line;
  if(++$ct == $step){
    CounterPrint("reads = $ct pairs = $pair_found");
    $step = $step + 100000;
  }
}
CounterPrint("\n");

sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}