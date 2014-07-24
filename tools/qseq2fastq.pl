#!/usr/bin/perl

use warnings;
use strict;

if($#ARGV<0){
   die "Usage: $0 <file>\n";
}

open(IN,$ARGV[0]) || die "Can't open $ARGV[0] for reading --fatal.\n";

while (<IN>) {
	chomp;
	my @parts = split /\t/;
	print "@";
        print "$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
	print "$parts[8]\n";
	print "+\n";
	print "$parts[9]\n";
}

close IN;
