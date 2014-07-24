  #############################################################
  #Marten Boetzer 13-06-2011                                  #
  #SSPACE perl subscript readLibFiles.pl                      #
  #This script;                                               #
  #  -reads, converts and filters original input sequences    #
  #############################################################

  use Storable;
  use File::Path;
  use File::Basename;
  use threads;

  my $seplines = ("-" x 60)."\n";
  my $maxlen = 0;

  my $libraryfile = $ARGV[0];
  my $base_name = $ARGV[1];
  my $extending = $ARGV[2];
  my $unpaired_file = $ARGV[3];
  my $min_overlap = $ARGV[4];
  my $thread = $ARGV[5];
  my $log = $base_name . ".logfile.txt";
  my $summaryfile = $base_name.".summaryfile.txt";

  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";

  my $filenameOutFilt = "filtered.readpairs.fasta";
  my $filenameOutExt = $base_name . ".singlereads.fasta";

#-------------------------------------------------READ UNPAIRED FILE CONTAINING SINGLE READS
  &readUnpairedFile($unpaired_file) if ($unpaired_file);
#-------------------------------------------------LOOP THROUGH EACH LIBRARY IN LIBRARYFILE AND STORE AND FILTER READS
  open(FILELIB, "< $libraryfile");

  my ($library, $fileA, $fileB, $insert_size, $insert_stdev, $reverse, $libResHash);
  my ($prevlibrary, $ctlib) = ("",0);
  &printMessage("\n=>".getDate().": Reading, filtering and converting input sequences of library file initiated\n");

  while(<FILELIB>){
    chomp;
    ($library, $fileA, $fileB, $insert_size, $insert_stdev, $reverse) = split(/\s+/, $_);

    next if($library eq "");
    $ctlib=0 if($library ne $prevlibrary && $prevlibrary ne "");
    $ctlib++;

    my ($fileBaseName1, $dirName1, $fileExtension1) = fileparse($fileA);
    my ($fileBaseName2, $dirName2, $fileExtension2) = fileparse($fileB);

    my $fname = "reads/$base_name.$library.filtered.readpairs.singles.fasta";
    my ($counter2, $Ncount2);
    #Process multiple files at the same time if multithreaded option is set (-T parameter larger than 1)
    if($fileA ne "TAB"  && $thread > 1){
       my $thr = threads->create(\&generateInputFiles, $library, $fileA, $fileB, $extending, $reverse, $fname, $ctlib);
       if(!($ctlib % $thread)){
         foreach my $thr (threads->list()) {
           my @res = $thr->join();
           ($lib,$nreads,$ncount) = split(/,/,$res[0]);
           $libResHash->{$lib}{'reads'}+=$nreads;
           $libResHash->{$lib}{'N'}+=$ncount;
         }
       }
    #otherwise, process only one file at a time
    }elsif($fileA ne "TAB" && $thread <=1){
      my $out = &generateInputFiles($library, $fileA, $fileB, $extending, $reverse, $fname, $ctlib);
      ($lib,$nreads,$ncount) = split(/,/,$out);
      $libResHash->{$lib}{'reads'}+=$nreads;
      $libResHash->{$lib}{'N'}+=$ncount;
    }
    #if user has inserted a TAB file, calculate read statistics
    if($fileA eq "TAB"){
      open FILE, "$fileB" or die $!;
      my ($fileBaseName2, $dirName2, $fileExtension2) = fileparse($fileB);
      print "Reading tabfile: $fileBaseName2...\n";
      $counter2++ while(<FILE>);
      $libResHash->{$lib}{'reads'}+=$counter2;
      $libResHash->{$lib}{'N'} = 0;
      close FILE;
    }
    $prevlibrary = $library;
  }
  #Process remaining reads
  if($fileA ne "TAB"){
    foreach my $thr (threads->list()) {
      my @res = $thr->join();
      ($lib,$nreads,$ncount) = split(/,/,$res[0]);
      $libResHash->{$lib}{'reads'}+=$nreads;
      $libResHash->{$lib}{'N'}+=$ncount;
    }
  }
  #Print read statistics to the summary file
  &printMessage("\n$seplines");
  foreach my $libs (keys %$libResHash){
    my $totcounter = $libResHash->{$libs}{'reads'};
    my $totNcount = $libResHash->{$libs}{'N'};
    my $filt = $totcounter-$totNcount;
    print SUMFILE "READING READS $libs:\n";
    print SUMFILE "$seplines\tTotal inserted pairs = $totcounter \n";
    print SUMFILE "\tNumber of pairs containing N's = $totNcount \n\tRemaining pairs = $filt\n$seplines\n";
  }
  close FILELIB;
  close SUMFILE;
  close LOG;

  mkpath('process_OK'); #make directory, indicating that process has run OK

#--------------------------------------------------

###CONVERT INPUT SEQUENCES BY REMOVING PAIRED READS HAVING AN 'N'
sub generateInputFiles{
  my ($lib, $fileA, $fileB, $extension, $reverse, $fname, $libct) = @_;
  my ($name,$seq1,$seq2, $res1,$res2);
  my ($counterext, $Ncount, $countsinglet, $fastq, $step) = (0,0,0,0,1000000);
  open (OUTSINGLEFILE, ">reads/$base_name.$lib.file$libct.fa") || die "Can't write to single file file$fname-- fatal\n";

  #check if file is fastQ or fastA
  open(TEST, "< $fileA");
  $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);

  open(FILEA, "< $fileA");
  open(FILEB, "< $fileB");
  CounterPrint("Reading read-pairs $lib.$libct @ $countsinglet       ");
  while(<FILEA>) {
    <FILEB>;
    $seq1 = uc(<FILEA>), $seq1 =~ s/^\r\n/\n/;
    $seq2 = uc(<FILEB>), $seq2 =~ s/^\r\n/\n/;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);

    $res1 = index($seq1,"N");
    $res2 = index($seq2,"N");
    #if both reads contain N's, do not use them for contig extension and for scaffolding
    if($res1 == -1 && $res2 == -1){
       print OUTSINGLEFILE ">read$countsinglet/1\n$seq1>read$countsinglet/2\n$seq2";
    }else{
      $Ncount++;
    }
    if(++$countsinglet == $step){   
      CounterPrint("Reading read-pairs $lib.$libct @ $countsinglet         ");
      $step = $step + 1000000;
    }

  }
  CounterPrint("\n") if($thread <= 1);
  CounterPrint((" " x 40));
  close OUTSINGLEFILE;
  close FILEB;
  close FILEA;
  return "$lib,$countsinglet,$Ncount";
}

#------------------READ UNPAIRED SINGLE READS FILE WHEN -u IS SET

sub readUnpairedFile{
  my ($file) = @_;
  open(INUNPAIRED, "< $file") || die "Can't open $file -- fatal\n";
  open OUTFILEExt, "> reads/$filenameOutExt";

  &printMessage("\n=>".getDate().": Reading, filtering and converting unpaired input sequences initiated\n");

  my ($seq1, $name);
  my ($counterext, $counter, $step, $fastq) = (0,0, 100000,0);

  open(TEST, "< $file");
  $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);
  while(<INUNPAIRED>) {
    $seq1 = uc(<INUNPAIRED>); $seq1 =~ s/\r\n/\n/; chomp $seq1;

    #FASTQ FORMAT
    if ($fastq){
      <INUNPAIRED>; <INUNPAIRED>;
    }
    # ELSE FASTA FORMAT
    if(index($seq1, "N") == -1){
       print OUTFILEExt ">$counterext\n$seq1\n";
       $counterext++;
    }
    if(++$counter == $step){
      CounterPrint($counter);
      $step = $step + 100000;
    }
  }
  CounterPrint("                ");

  print SUMFILE "READING UNPAIRED READS:\n";
  print SUMFILE "$seplines\tTotal inserted reads = $counter \n";
  print SUMFILE "\tNumber of reads containing N's = ".($counter-$counterext)."\n\tRemaining reads = $counterext\n";
  close OUTFILEext;
  close INUNPAIRED;
}

###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGC/TACG/;
   return (reverse());
}

###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
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

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(SUMFILE), $| = 1)[0]);
  select((select(LOG), $| = 1)[0]);
  $|++;
}

#########END readLibFiles.pl