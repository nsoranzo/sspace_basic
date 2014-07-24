#!/usr/bin/perl

use strict;

if($#ARGV<0){
   die "Usage: $0 <file>\n";
}

open(IN,$ARGV[0]) || die "Can't open $ARGV[0] for reading --fatal.\n";
my $fasta = $ARGV[0] . ".fa";
open(OUT,">$fasta") || die "Can't open $fasta for writing --fatal.\n";

while (<IN>) {
	chomp;
	my @parts = split(/\s+/);
        my $concat = ">$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]";
        print OUT "$concat\n";
	print OUT "$parts[8]\n";
}

close IN;
close OUT;

exit;
