#!/usr/bin/perl

#AUTHOR
# Marten Boetzer and Walter Pirovano (c) 2010
# SSAKE-based Scaffolding of Pre-Assembled Contigs after Extension (SSPACE)
# walter.pirovano@baseclear.com

#NAME
#   SSPACE Marten Boetzer - Walter Pirovano November 2011

#SYNOPSIS
#   SSAKE-based Scaffolding of Pre-Assembled Contigs after Extension (SSPACE)

#DOCUMENTATION
#   README, MANUAL and TUTORIAL distributed with this software @ www.baseclear.com
#   Boetzer M, Henkel VJ, Jansen HJ, Butler D and Pirovano W. 2011. Scaffolding pre-assembled contigs using SSPACE. Bioinformatics 27(4) p578-9.
#   http://www.baseclear.com/sequencing/data-analysis/bioinformatics-tools/
#   We hope this code is useful to you -- Please send comments & suggestions to Walter.Pirovano@baseclear.com
#   If you use either the SSPACE code or ideas, please cite our work appropriately and accurately

#LICENSE
#   SSPACE Copyright (c) 2010-2011 BaseClear B.V. The Netherlands. All rights reserved.
#   SSAKE Copyright (c) 2006-2010 Canada's Michael Smith Genome Science Centre. All rights reserved.

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   note: insert size and distance between pairing reads are used interchangeably

#MAJOR CHANGES ON SSAKE V3.4 TO FORM SSPACE
#   -New scaffolding feature dealing with contigs having multiple alternatives
#   -Seperate scripts to decrease the memory usage
#   -Automatic filtering of reads and duplicate mate pairs
#   -Option for contig extension on unfiltered reads
#           -Removed tracking of reads during contig extension
#           -Mapping single reads to extended and non extended contigs
#   -Single reads mapped more than once to a contig are removed for scaffolding
#   -A summary file is generated containing detailed information about the scaffolding process
#   -An evidence file is generated which indicates the contigs present in the scaffolds
#   -Optional; Scaffolds and their contigs are visualised by generating a .dot file

#MAJOR CHANGES ON SSPACE Basic v2.0;

# GENERAL
#   -Last column of the library file should be the orientation of the reads, instead of indication of being reverse complement or not. Options are FR, FF, RF and RR.
#   -Fixed some bugs in the summary file and removed some useless information.
#   -Included the -z option which specifies the minimal contig length that will be used for scaffolding. Contigs below this length are discarded for scaffolding.
#   -Included the possibility to include TAB delimited files with read mapping information, format is; ctg1 start1 end1 ctg2 start2 end2
#             - if a read is reverse complement on a contig, start and end should be turned around e.g. ctg1 100 150 ctg2 150 100 indicates that the second read is reverse complement on ctg2
#             - No contig filtering can be applied if TAB delimited files are included
#             - See MANUAL for more information of how to use the tab file option
#   -Included some scripts to convert a .sam file to a .tab file

# BOWTIE
#   -Included the -g option to specify maximum allowed gaps for Bowtie. This option corresponds to the -v option in Bowtie.
#   -Now able to do multithreaded Bowtie using the -T option (-T 3 does 3 threads). This option corresponds to the -p option in Bowtie.

# READING FILES:
#   -Speeded up the reading of the library files for a single threaded run
#   -Now able to read multiple libraries at once using the multithread -T option. -T 3 reads three files at the same time.

# CONTIG EXTENSION
#   -Included the -r option for contig extension (default is 0.9).
#   -Speeded up and reduced the memory usage during the contig extension.
#             - SSPACE reads in the output of Bowtie at once, rather than reading it from the output file.
#             - Faster check for presence of subsequence of a read, thereby able to faster check for overlapping sequences with the contig.

# SCAFFOLDING
#   -Combined the functions readBowtie and pairContigs, which saves runtime and memory.
#   -Saving runtime by reading Bowtie results in at once, instead of reading it from Bowtie's output file.
#   -Included a pre-filtering step of multiple alternative contig links before scaffolding. This step was previously done during scaffolding, now it's a step before scaffolding. It reduces the number of errors within the scaffolds.
#   -Additional check to connect two alternative contigs, making the scaffolds more reliable, especially with mate pair libraries. The search space is included in the calculation of the ratio, rather than looking at the number of links only. See the README file for more information.
#   -Calculation of mean insert size based on mapped read pairs on same contig. Users can choose this value for better estimation of gap sizes. Especially for paired-end sequences.

#   -Fixed a bug in the mergeContigs function. Indication of contigs merged in previous libraries were not displayed in the final .evidence file.

#-------------------------------------------------LOAD PACKAGES AND DEFINE VARIABLES
  use strict;
  use Storable;
  use Getopt::Std;
  use File::Path;
  use File::Basename;

  #Specify path to DotLib
  use FindBin qw($Bin);
  use lib "$Bin/dotlib/";
  use DotLib;

  use vars qw($opt_m $opt_o $opt_v $opt_p $opt_k $opt_a $opt_z $opt_s $opt_b $opt_n $opt_l $opt_x $opt_u $opt_t $opt_T $opt_g $opt_r);
  getopts('m:o:v:p:k:a:z:s:b:n:l:x:u:t:T:g:r:');
  my ($base_overlap, $min_overlap, $verbose, $MIN_READ_LENGTH, $SEQ_SLIDE, $min_base_ratio, $min_links, $max_link_ratio, $unpaired_file, $max_trim, $base_name, $max_count_trim, $min_tig_overlap, $doplot, $extending, $threads, $minContigLength, $gaps, $unpaired, $gapclosure) = (20, 32, 0, 16, 1, 0.9, 5, 0.7, "no-u", 0, "standard_output", 10, 15, 0, 0, 1, 0, 0, 0, 0);

  my $version = "[SSPACE_Basic v2.1]";
  my $seplines = ("-" x 60)."\n";
  my ($MAX, $MAX_TOP, $TRACK_COUNT) = (0, 100, 1);# $MAX_TOP is the very maximum anchoring edge sequence that will be searched

#-------------------------------------------------READ OPTIONS

  if(!($opt_l) || !($opt_s)){
    print STDERR "ERROR: Parameter -l is required. Please insert a library file\n" if(!$opt_l);
    print STDERR "ERROR: Parameter -s is required. Please insert a contig FASTA file\n" if(!$opt_s);
    my $error_msg = <<"END_MSG";
\nUsage:\n
============ General Parameters ============\n
-l  Library file containing two paired read files with insert size, error and either mate pair or paired end indication (REQUIRED)\n
-s  FASTA file containing contig sequences used for extension. Inserted pairs are mapped to extended and non-extended contigs (REQUIRED)\n
-x  Indicate whether to extend the contigs of -s using paired reads in -l (-x 1=extension, -x 0=no extension, default -x $extending)\n
============ Extension Parameters ============\n
-m  Minimum number of overlapping bases with the seed/contig during overhang consensus build up (default -m $min_overlap)\n
-o  Minimum number of reads needed to call a base during an extension (default -o $base_overlap)\n
-t  Trim up to -t base(s) on the contig end when all possibilities have been exhausted for an extension (default -t $max_trim)\n
-u  FASTA/FASTQ file containing unpaired sequence reads (optional)\n
-r  Minimum base ratio used to accept a overhang consensus base (default -r $min_base_ratio)\n
============ Scaffolding Parameters ============\n
-z  Minimum contig length used for scaffolding. Filters out contigs below this value (default -z $minContigLength)\n
-k  Minimum number of links (read pairs) to compute scaffold (default -k $min_links)\n
-a  Maximum link ratio between two best contig pairs. *Higher values lead to least accurate scaffolding* (default -a $max_link_ratio)\n
-n  Minimum overlap required between contigs to merge adjacent contigs in a scaffold (default -n $min_tig_overlap)\n
============ Bowtie Parameters ============\n
-g  Maximum number of allowed gaps during mapping with Bowtie. Corresponds to the -v option in Bowtie. *Higher number of allowed gaps can lead to least accurate scaffolding* (default -g $gaps)\n
-T  Specify the number of threads in Bowtie. Corresponds to the -p/--threads option in Bowtie (default -T $threads)\n
============ Additional Parameters ============\n
-b  Base name for your output files (default -b $base_name)\n
-v  Runs in verbose mode (-v 1=yes, -v 0=no, default -v $verbose)\n
-p  Make .dot file for visualisation (-p 1=yes, -p 0=no, default -p $doplot)
END_MSG
    die $error_msg;
  }

  my $libraryfile = $opt_l if ($opt_l);
  my $filecontig = $opt_s if($opt_s);
  $extending = $opt_x if($opt_x eq 1);
  $min_overlap = $opt_m if ($opt_m);
  $base_overlap = $opt_o if ($opt_o);
  $max_trim = $opt_t if ($opt_t);
  $unpaired_file = $opt_u if($opt_u);
  $min_base_ratio = $opt_r if ($opt_r);
  $minContigLength = $opt_z if($opt_z);
  $min_links = $opt_k if ($opt_k);
  $max_link_ratio = $opt_a if ($opt_a);
  $min_tig_overlap = $opt_n if($opt_n);
  $gaps = $opt_g if($opt_g);
  $threads = $opt_T if ($opt_T);
  $base_name = $opt_b if($opt_b);
  $verbose = $opt_v if ($opt_v);
  $doplot = $opt_p if($opt_p);

#-------------------------------------------------CHECKING PARAMETERS
  die "ERROR: Invalid (-l) library file $libraryfile ...Exiting.\n" if(! -e $libraryfile);
  die "ERROR: Invalid (-s) contig file $filecontig ...Exiting.\n" if(! -e $filecontig);
  die "ERROR: -x must be either 0 or 1. Your inserted -x is $extending...Exiting.\n" if(!($extending == 0 || $extending == 1));
  die "ERROR: -m must be a number between 15-50. Your inserted -m is $min_overlap ...Exiting.\n" if(!($min_overlap =~ /^\d+$/) || $min_overlap < 10 || $min_overlap > 50);
  die "ERROR: -o must be set to 1 or higher. Your inserted -o is $base_overlap ...Exiting.\n" if($base_overlap < 1);
  die "ERROR: -t must be a positive integer. Your inserted -t is $max_trim ...Exiting.\n" if(!($max_trim =~ /^\d+$/));
  die "ERROR: Invalid unpaired file $unpaired_file -- fatal\n" if(! -e $unpaired_file && $opt_u);
  die "ERROR: -r must be a number between 0.0 and 1.0. Your inserted -r is $min_base_ratio ...Exiting.\n" if($min_base_ratio < 0 || $min_base_ratio > 1);
  die "ERROR: -z must be a positive integer. Your inserted -z is $minContigLength...Exiting.\n" if (!($minContigLength =~ /^\d+$/));
  die "ERROR: -k must be a positive integer. Your inserted -k is $min_links ...Exiting.\n" if(!($min_links =~ /^\d+$/));
  die "ERROR: -a must be a number between 0.0 and 1.0. Your inserted -a is $max_link_ratio ...Exiting.\n" if($max_link_ratio < 0 || $max_link_ratio > 1);
  die "ERROR: -n must be a positive integer. Your inserted -n is $min_tig_overlap ...Exiting.\n" if (!($min_tig_overlap =~ /^\d+$/));
  die "ERROR: -g must be a positive integer between 0 and 3. Your inserted -g is $gaps...Exiting.\n" if (!($gaps =~ /^\d+$/) || $gaps > 3);
  die "ERROR: -T must be a positive integer. Your inserted -T is $threads...Exiting.\n" if (!($threads =~ /^\d+$/));
  die "ERROR: -p must be either 0 or 1. Your inserted -p is $doplot...Exiting.\n" if(!($doplot == 0 || $doplot == 1));

#-------------------------------------------------check library file;
  open(FILELIB, "< $libraryfile");
  my ($min_allowed, $library, $fileA, $fileB, $insert_size, $insert_stdev, $orientation);
  my $countline=0;
  while(<FILELIB>){
    chomp;
    $countline++;
    my @line = split(/\s+/, $_);
    if($#line >= 0){
      if($opt_l){
        die "ERROR: Line $countline in your library file ($libraryfile) contains $#line spaces, which should be 5 spaces. Check that no spaces are within the file names.\n" if($#line != 5);

        my ($library, $fileA, $fileB, $insert_size, $insert_stdev, $orientation) = split(/\s+/, $_);
        if($fileA ne "TAB"){
          die "ERROR: Invalid file in library $library: $fileA -- fatal\n" if(! -e $fileA);
        }else{
          die "ERROR: Can't apply filtering using the -z option (-z = $minContigLength)and insertion of a TAB file -- fatal\n" if($minContigLength > 0);
        }
        die "ERROR: Invalid file in library $library: $fileB -- fatal\n" if(! -e $fileB);
        die "ERROR: Insert size should be higher than or equal to 0. Your library $library has insert size of $insert_size. Exiting.\n" if(!($insert_size>0) || !($insert_size =~ /^\d+$/));
        die "ERROR: Insert stdev must be a number between 0.00 and 1.00. Your library $library has insert size of $insert_stdev. Exiting.\n" if($insert_stdev < 0 || $insert_stdev > 1 || !($insert_stdev * 1 eq $insert_stdev));
        die "ERROR: Orientation must have length of 2 characters and should contain one of the following; FR, FF, FR or RF. Your library $library has orientation of $orientation ...Exiting.\n" if(!(length($orientation) == 2) || !($orientation =~ /[FR][FR]/));
      }
    }
  }
  close FILELIB;
#-------------------------------------------------Make folder structure
  mkpath('intermediate_results');
  mkpath('pairinfo');
  mkpath('reads');
  mkpath('bowtieoutput');

  $unpaired = $unpaired_file if (-e $opt_u && $extending == 1);
#-------------------------------------------------Print input parameters
  my $contig = "intermediate_results/" . $base_name .  ".formattedcontigs.fasta";

  my $log = $base_name . ".logfile.txt";
  my $summaryfile = $base_name.".summaryfile.txt";
  open (LOG, ">$log") || die "Can't write to $log -- fatal\n";

  open (SUMFILE, ">$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  close SUMFILE;

  my $init_message =  "Your inserted inputs on $version at ".getDate().":\nRequired inputs: \n\t-l = $libraryfile\n\t-s = $filecontig\n\t-b = $base_name\n\n";
  $init_message .= "Optional inputs:\n\t-x = $extending\n\t-z = $minContigLength\n\t-k = $min_links\n";
  $init_message .=  "\t-a = $max_link_ratio\n\t-n = $min_tig_overlap\n\t-T = $threads\n\t-p = $doplot\n\n";

  $init_message .= "Contig extension inputs:\n\t-o = $base_overlap\n\t-t = $max_trim\n\t-m = $min_overlap\n\t-r = $min_base_ratio\n\n" if($extending == 1);

  &printMessage($init_message);
  close LOG;
#-------------------------------------------------READING AND CONVERTING INPUT SEQUENCES
  system("perl $Bin/bin/readLibFiles.pl $libraryfile $base_name $extending $unpaired $min_overlap $threads");
  checkStatus();
#-------------------------------------------------FORMATTING OR EXTENDING CONTIGS
  system("perl $Bin/bin/ExtendOrFormatContigs.pl $contig $base_name $extending $filecontig $MIN_READ_LENGTH $base_overlap $min_overlap $min_base_ratio $max_trim $verbose $Bin $minContigLength $libraryfile $gaps $threads");
  checkStatus();
#--------------------------------------------------UPDATE SUMMARY FILE
  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";

  #write summary of initial contigs
  my $sumfile .= "\nSUMMARY: \n".$seplines."\tInserted contig file;\n";
  $sumfile = &writesummaryfiles($filecontig, "contig", $sumfile);
  #write summary of extended contigs
  my $extended_tig = "intermediate_results/" . $base_name .  ".extendedcontigs.fasta";
  $sumfile .= "\tAfter extension;\n" if($extending);
  $sumfile = &writesummaryfiles($extended_tig, "contig", $sumfile) if($extending);

  #write summary of filtered contigs
  if($minContigLength > 0){
    $sumfile .= "\tAfter filtering (z >= $minContigLength);\n";
    $sumfile = &writesummaryfiles($contig, "contig", $sumfile);
  }else{
    $contig = $extended_tig if($extending);
  }
  &FlushFiles();
  close LOG;
  close SUMFILE;

#--------------------------------------------------GO THROUGH EACH LIBRARY AND SCAFFOLD
  open(FILELIB, "< $libraryfile") || die "Can't open $libraryfile -- fatal\n";
  my ($lib, $fileA, $fileB, $insert_size, $insert_stdev, $pair, $headscaffolds, $prevlib, $mergedtigs, $evidencefile);

  while(<FILELIB>){
    chomp;
    &FlushFiles();
    ($lib, $fileA, $fileB, $insert_size, $insert_stdev, $orientation) = split(/\s+/, $_);
    next if($lib eq $prevlib || $lib eq '');

    my $tabfile = 0;
    $tabfile = 1 if($fileA eq "TAB");

    $prevlib = $lib;
    $min_allowed = -1 * ($insert_stdev * $insert_size);

    open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
    &printMessage("\nLIBRARY $lib\n".$seplines);
    close LOG;

    open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
    print SUMFILE "\n\nLIBRARY $lib STATS:\n".("#" x 80),"\n";
    close SUMFILE;

    my $scaffold = "intermediate_results/" . $base_name . ".$lib.scaffolds";
    $mergedtigs = "intermediate_results/" . $base_name . ".$lib.scaffolds.fasta";
    my $issues = "pairinfo/" . $base_name . ".$lib.pairing_issues";
    my $distribution = "pairinfo/" . $base_name . ".$lib.pairing_distribution.csv";

#-------------------------------------------------MAPPING READ PAIRS USING FILTERED FASTA FILE
    mkpath("tmp.$base_name");
#-------------------------------------------------Scaffold the contigs and generate .scaffold file
    system("perl $Bin/bin/PairingAndScaffolding.pl $Bin $gaps $contig $base_name $issues $distribution $verbose $lib $insert_size $min_allowed $scaffold $min_links $max_link_ratio $orientation $threads") if(!$tabfile);
    system("perl $Bin/bin/PairingAndScaffolding.pl $Bin $gaps $contig $base_name $issues $distribution $verbose $lib $insert_size $min_allowed $scaffold $min_links $max_link_ratio $orientation $threads $tabfile $fileB $filecontig $evidencefile") if($tabfile);
    checkStatus();

    #retrieve the contigs that were stored
    my $contigstored = "tmp.$base_name/contigs.stored";
    my $contigs = retrieve("$contigstored");
#-------------------------------------------------Generate .fasta file and .evidence file with scaffolds
    open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
    ($headscaffolds, $evidencefile) = &mergeContigs($scaffold, $contigs, $mergedtigs, 50, $verbose, $min_tig_overlap,$max_count_trim);
    $contig = $mergedtigs;
#-------------------------------------------------write summary of scaffolds
    $sumfile .= "\tAfter scaffolding $lib:\n";
    $sumfile = &writesummaryfiles($mergedtigs, "scaffold", $sumfile);

#-------------------------------------------------
    open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
    print SUMFILE ("#" x 80),"\n";
    close SUMFILE;
    &printMessage("\n$seplines");
    $contigs = (''); undef $contigs;

    my $removedir = "tmp.$base_name";
    rmtree([$removedir, 'blurfl/quux']);  #remove 'tmp' folder
  }#END OF LIBRARY LOOP

  #-------------------------------------------------END OF LIBRARIES. PRINT SUMMARY TO FILE AND END SESSION
  my $finalfile = $base_name . ".final.scaffolds.fasta";
  my $finalevfile = $base_name . ".final.evidence";

  open (EVID, $evidencefile);
  open (FINALEV, "> $finalevfile");
  while(<EVID>){
    print FINALEV $_;
  }

  open (SCAF, $mergedtigs);
  open (FINAL, "> $finalfile");
  while(<SCAF>){
    print FINAL $_;
  }

  #make .dot file for visualisation
  &visualiseScaffolds($base_name.".visual_scaffolds", $evidencefile) if($doplot);

  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  &printMessage("\n=>".getDate().": Creating summary file\n");
  print SUMFILE $sumfile.$seplines;
  my $time = (time - $^T);
  my $minutes = int ($time / 60);
  $time = $time % 60;
  &printMessage(("*" x 50)."\n\nProcess run succesfully on ".getDate()." in $minutes"." minutes and $time"." seconds\n\n\n");
  close SCAF;
  close FINAL;
  close EVID;
  close FINALEV;
  close LOG;
  close SUMFILE;
  #END OF MAIN PROGRAM

###MAKE A .FASTA FILE OF THE FOUND SCAFFOLDS. EITHER MERGE TWO CONTIGS WHEN A OVERLAP OF -n EXISTS OR PLACE A GAP
sub mergeContigs{

   my ($scaffold, $contigs, $mergedtigs, $chunk, $verbose,$min_tig_overlap,$max_count_trim) = @_;

   &printMessage("\n=>".getDate().": Merging contigs and creating FASTA file of scaffolds\n");

   open(IN,$scaffold) || die "can't read $scaffold -- fatal\n";

   my $evidence_file = $mergedtigs;
   $evidence_file =~ s/.fasta/.evidence/;
   open(SCAFS,">$evidence_file") || die "can't write to $evidence_file -- fatal\n";
   open(OUT,">$mergedtigs") || die "can't write to $mergedtigs -- fatal\n";
   my $scafhashcount = keys ( %$headscaffolds );
   my $scaffoldHashStart;
   my ($tot,$sct,$ct_merge, $step) = (0,0,0,100);
   while(<IN>){### each line is a scaffold
      chomp;
      my @a = split(/\,/);
      my @tig;

      if($a[2]=~/\_/){
         @tig = split(/\_/,$a[2]);
      }else{
         push @tig, $a[2];
      }
      if(++$sct == $step){
        CounterPrint($sct);
        $step = $step + 100;
      }
      my ($ct,$tigsum,$mct,$prev,$word,$template,$seq,$prevseq,$headconcat,$prevEstimatedDistance, $prevLinks) = (0,0,0,"NA","NA","NA","","","","");
      foreach my $t (@tig){### each contig
         $ct++;

         if($t=~/([fr])(\d+)z(\d+)(\S+)?/i){

            my $orient = $1;
            my $tnum=$2;
            my $head = $orient . $tnum;
            my $search = "tig" . $tnum;
            my $other = $4;
            $tot+= $3;
            $tigsum +=$3;

            my ($estimatedDistance, $links) = ("", "");
            $estimatedDistance = $1 if($other=~/m((\-)?\d+)/);
            $links = $1 if($other=~/k((\-)?\d+)/);
            print "\tSC $a[0] - TIG $ct.  pattern: $t search: $search totalTigSize: $tot Orientation: $orient Gap/Overlap estimated distance: $estimatedDistance\n" if($verbose);

            my $count_trim = 0;

            $seq = $contigs->{$tnum}{'seq'};
            $seq = reverseComplement($seq) if($orient eq "r");
            chomp $seq;
            my $prev;
            if($scafhashcount >0){
              $prev = $headscaffolds->{$tnum}{'head'};
              $prev =~ s/^\n//;
              chomp $prev;
              delete $headscaffolds->{$tnum};
              chomp $prev;
              if($orient eq "r"){  ###Reverse all contigs if the whole scaffold is a reverse complement. ftig -> rtig and rtig -> ftig
                my @prevarray = split("\n", $prev);
                if($#prevarray >=0){
                  my $newprev="";
                  my ($tnum, $sizetig, $links, $gap, $prevline, $merge) = ("","","","","","");
                  for(my $i = $#prevarray; $i >= 0; $i--){

                    my @info = split(/\|/, $prevarray[$i]);
                    if($#info eq 1){
                      ($tnum, $sizetig) = split(/\|/, $prevarray[$i]);
                    }else{
                      ($tnum, $sizetig, $links, $gap, $merge) = split(/\|/, $prevarray[$i]);
                    }
                    $tnum =~ tr/fr/rf/;
                    if($prevline ne ""){
                      $newprev .= $prevline."|".$links."|".$gap."\n" if($merge eq "");
                      $newprev .= $prevline."|".$links."|".$gap."|".$merge."\n" if($merge ne "");
                    }
                   $prevline = $tnum."|".$sizetig;
                  }
                  $newprev .= $prevline;
                  $prev = $newprev;
                }
              }
            }
            else{
              $prev = "$orient"."_$search|size".length($seq);
            }
              $prev .= "|links$links|gaps$estimatedDistance" if($links ne "");


            #print "$prev\n";
            if($word ne "NA"){
               #####
               if(length($seq)<=$chunk){
                  $template = $seq;
               }else{
                  $template = substr($seq,0,$chunk);
               }

               ##### word search
               my $dynamic_word = $word;
	       if($prevEstimatedDistance <= 0){
                 SCAN:
                 until($template =~ /$dynamic_word/){
                   $dynamic_word = substr($dynamic_word,1,length($dynamic_word));
                   if(length($dynamic_word) < $min_tig_overlap){
                     $count_trim++;
                     last SCAN if($count_trim >= $max_count_trim);
                     $dynamic_word = substr($word,0,length($word)-$count_trim);
                   }
                 }
	       }
               if($prevEstimatedDistance <= 0  && $seq =~ /^\S{0,$max_count_trim}$dynamic_word(.*)/){### will grab the left-most match which is ok
                  my $tail = $1;
                  my $all = "ERROR_";
                  while($prevseq =~ /^(.*)$dynamic_word/ig){
                     $all = $1;
                  }
                  print "$prevseq **** $all **** WORD:$word *** DWord:$dynamic_word *** COUNTTRIM:$count_trim\n" if($all=~/ERROR/);

                  $prevseq = $all . lc($dynamic_word) . $tail;
                  my $overlap = length($dynamic_word);
                  $ct_merge++;
                  print "$ct_merge. GROUNDS FOR MERGING ($overlap nt overlap) !!!\n" if($verbose);
                  $headconcat .= "|merged$overlap"."\n".$prev;
               }else{
                  ### ADDED RLW 5.MAR.2010
                  if($prevEstimatedDistance <= 0){
                     $prevseq .= "n" . $seq
                  }else{
                     $prevseq .= ("N" x $prevEstimatedDistance) . $seq;
                  }
                  $headconcat .= "\n".$prev;

               }
            }else{
               $prevseq = $seq;
               $headconcat = "\n".$prev;
               $mct++;
            }

            ##### For the next search
            if(length($seq)<=$chunk){
               $word = $seq;
            }else{
               $word = substr($seq,length($seq)-$chunk,$chunk); ### this will be the next word to search with
            }
            $prevEstimatedDistance = $estimatedDistance;
            $prevLinks = $links;
         }#tig regex

      }#each tig
      my $scsz = length($prevseq);
      $scaffoldHashStart->{$sct}{'head'} = $headconcat;

      my @line = split(/\n/, $headconcat);
      print SCAFS ">$a[0]|size$scsz|tigs".($#line)."$headconcat\n\n";
      print OUT ">$a[0]|size$scsz\n$prevseq\n";
      $prevseq = '';
   }
   close IN;
   close SCAFS;
   close OUT;
   CounterPrint("                ");
   undef $contigs;
   &FlushFiles();
   return ($scaffoldHashStart, $evidence_file);
}
###WRITE SUMMARY STATISTICS FOR ALL CONTIGS OR SCAFFOLDS
sub writesummaryfiles{
  my ($input_file, $insert, $sumfile) = @_;

  open (INFILE, $input_file) || die "Can't open input file $input_file.\n";

  my ($counter, $sum, $seq, $name, $foundN50, $sumN50, $totalNcount) = (0,0, "","", 0, 0);
  my (@line, @lengths);
  while (<INFILE>) {
    s/\r\n/\n/;
    chomp;
    $seq.= $_ if(eof(INFILE));
    if ($_ =~ /^[>]/ || eof(INFILE)) {
      if($counter > 0){
         push(@lengths, length($seq));
         $sum+= length($seq);
         my $Ncount = () = $seq =~ /[Nn]/g;
         $totalNcount += $Ncount;
         ($seq) = "";
      }
      $counter++;
    }
    else {
       $seq .= $_;
    }
  }
  $counter--;
  my $half_length = $sum/2;

  my @lengths2 = reverse sort { $a <=> $b } @lengths;

  for(my $i = 0; $i <= $#lengths && $foundN50 == 0; $i++)
  {
    $sumN50 += @lengths2[$i];
    if($sumN50 >= $half_length){
      $foundN50 = @lengths2[$i] if($sumN50 >= $half_length);
      last;
    }
  }
  $sumfile .= "\t\tTotal number of $insert"."s = $counter\n";
  $sumfile .= "\t\tSum (bp) = ". $sum. "\n";
  $sumfile .= "\t\t\tTotal number of N's = $totalNcount\n";
  $sumfile .= "\t\t\tSum (bp) no N's = ". ($sum-$totalNcount)."\n";
  $sumfile .= "\t\tMax $insert size = ". @lengths2[0]."\n";
  $sumfile .= "\t\tMin $insert size = ". @lengths2[$#lengths]."\n";
  $sumfile .= "\t\tAverage $insert size = ".int($sum/$counter)."\n";
  $sumfile .= "\t\tN50 = ". $foundN50. "\n\n";

  close (INFILE);
  close OUTFILE;

  return $sumfile;
}


###FUNCTION TO GENERATE A VISUALISATION OF THE SCAFFOLDS AND THEIR CONTIGS IN .DOT FORMAT
sub visualiseScaffolds{
   my ($dotname, $evidence) = @_;
   my ($filext, $sizecutoff) = (1, 5000000);
   mkpath('dotfiles');
   my $filename2 = "dotfiles/$dotname.part".$filext.".dot";
   &printMessage("\n=>".getDate().": Producing .dot file for visualisation\n");

   open(IN,$evidence) || die "can't read $evidence -- fatal\n";
   open(DOT, ">$filename2") || die "can't open $filename2 -- fatal\n";
   printHeader(\*DOT, undef);
   my ($prevtig, $prevgap, $prevlinks, $prevratio, $scafcount) = ("","","", "",0);
   while(<IN>){
     chomp;
     my $line = $_;
     my $filesize = -s $filename2;

     if ($line =~ /^[>]/){
      endCluster(\*DOT) if($scafcount > 0);
       my $filesize = -s $filename2;
       if($filesize > $sizecutoff){
         printFooter(\*DOT);
         close(DOT);
         $filext++;
         $filename2 = "$dotname.part".$filext.".dot";
         open(DOT, ">$filename2") || die "can't open $filename2 -- fatal\n";
         printHeader(\*DOT, undef);
       }
       $scafcount++;
       $line =~ tr/[>\|]/ /;
       startCluster(\*DOT, $scafcount, "$line");
       ($prevtig, $prevgap, $prevlinks, $prevratio) = ("","","", "");
     }
     elsif($line =~ /^[fr]/){
        my @info = split(/\|/, $line);
        my ($tnum, $sizetig, $links, $gap);
        if($#info eq 1){
          ($tnum, $sizetig) = split(/\|/, $line);
        }else{
          ($tnum, $sizetig, $links, $gap) = split(/\|/, $line);
        }
        my ($orient, $tig) = split(/_/,$tnum);
        my $ori=-1;
        my ($other, $gap2) = split(/gaps/,$gap);
        my ($other, $links2) = split(/links/,$links);
        $ori = 1 if($orient eq "f");
        printNode(\*DOT, $tig, "$tig ($sizetig)", $ori);
        printEdge(\*DOT, $prevtig, $tig, "gap = $prevgap links = $prevlinks", undef) if($prevtig ne "");

        $prevtig = $tig;
        $prevgap = $gap2;
        $prevlinks = $links2;
     }
   }
   endCluster(\*DOT) if($scafcount > 0);
   printFooter(\*DOT);
   close(DOT);
   close IN;
}


###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGC/TACG/;
   return (reverse());
}

###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
  print LOG $message;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}

###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(SUMFILE), $| = 1)[0]);
  select((select(LOG), $| = 1)[0]);
  $|++;
}
#########END MAIN SCRIPT


sub checkStatus{
  &printMessage(("*" x 50)."\n\nProcess failed on ".getDate()."\n\n\n"), exit if(!(-d "process_OK"));
  rmtree(["process_OK", 'blurfl/quux']);
}
