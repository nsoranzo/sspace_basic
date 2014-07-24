#!/usr/bin/env python

__doc__ = """
TQS

Trim Quality Solexa Sequences (TQS)

SYNOPSIS
   Quality trim solexa-Illumina sequence reads using user-defined thresholds 
"""
__author__ = "Rene L. Warren"
__version__ = '1.0'

#LICENSE
#   Copyright (c) 2007 Canada's Michael Smith Genome Science Centre.  All rights reserved.

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

import sys, os, re, string
from datetime import datetime
from optparse import OptionParser


def main():
	usage = "Usage: %s --help"

	parser = OptionParser()
	parser.add_option("-f", "--sequence file", dest="seqfile",
	                  help="Illumina sequence file - Output format from the 1G Genome Analyzer (_seq.txt):                                       7       1       255     669     AACCCCCACTCCTACAACGCCATCATTCCCCTCGAC",)
	parser.add_option("-q", "--qual file", dest="qualfile",
	                  help="A prb file containing all the Illumina intensities, as outputted by the 1G Genome Analyzer (_prb.txt)",)
	parser.add_option("-l", "--length", dest="mer", type="int", default=36,
	                  help="Length of sequence reads (i.e. Number of sequencing cycles, default=36)",)
        parser.add_option("-t", "--threshold", dest="threshold", type="int", default=5,
                          help="Base intensity threshold value (-40 to 40, default=5)",)
        parser.add_option("-d", "--difference", dest="diff", type="int", default=5,
                          help="Base intensity difference between top intensity and second best (1 to 80, default=5)",)
        parser.add_option("-c", "--consec", dest="consec", type="int", default=20,
                          help="Minimum number of consecutive bases passing threshold values (default=20)",)
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
	                  help="Runs in Verbose mode.",)
	(opts, args) = parser.parse_args()
	
	try:
		f = open(opts.seqfile)
		seq = f.readlines()
		f.close()
	except Exception, e:
		print "ERROR: Could not read from %s: %s" % (opts.seqfile, e)
		print usage % (sys.argv[0:])
		sys.exit()

        try:
                f = open(opts.qualfile)
                qual = f.readlines()
                f.close()
        except Exception, e:
                print "ERROR: Could not read from %s: %s" % (opts.qualfile, e)
                print usage % (sys.argv[0:])
                sys.exit()
	

	fasta = "%s_I%sD%sL%s.trim.fa" % (opts.seqfile,opts.threshold,opts.diff,opts.consec)
	log = "%s.log" % opts.seqfile


        try:
                FASTA = open(fasta, 'w')
        except:
                print "ERROR: Can not write to %s" % fasta
                sys.exit()

	try:
		LOG = open(log, 'w')
	except:
		print "ERROR: Can not write to %s" % log
		sys.exit()
	

	if opts.mer < 15 or opts.mer > 200:
		print "ERROR: -l must be a number between 15 and 200."
		sys.exit()
	
	if opts.consec < 16 or opts.consec > opts.mer:
		print "ERROR: -c must be a number between 16 and -l."
		sys.exit()

	LOG.write("""
Running:
%s
-f %s
-q %s
-l %s
-c %s
-t %s
-d %s
Fasta file: %s

""" % (sys.argv[0:],opts.seqfile, opts.qualfile, opts.mer, opts.consec, opts.threshold, opts.diff, fasta))
	
	t0 = datetime.now()
	LOG.write("\nReading Quality File: %s\n" % str(t0)[:len('2006-10-05 23:04')])
	trim_info = parseQualFile(opts.threshold, opts.diff, opts.consec, opts.mer, qual, opts.verbose, LOG)
        t1 = datetime.now()
        LOG.write("\n\nTrimming low quality bases: %s\n" % str(t1)[:len('2006-10-05 23:04')])
	readNTrim(trim_info, seq, opts.verbose, FASTA, LOG)
        LOG.write("DNA sequences have been trimmed accordingly and placed in %s" % fasta)
	
	LOG.close()
	FASTA.close()
	return	

#--------------------------------------------------------------------------------------
def parseQualFile(threshold, difference, consecutive, read_length, qual, verbose, LOG):
	"""
	Parse a solexa-illumina intensity file
	
	Return a Dictionary of sequence order number, with the index value and length to extract 
	"""
	trim_info = {}
	ok_read = 0
	read_number = 0

	if verbose:
		print "Printing trimming pattern for all reads passing the set threshold values...\n"
	
	for line in qual:
		read_number += 1                ### this keeps track of the read order, respected between the prb and seq files
		concat = ""			### concat builds a string of bases passing the user-defined filter 
		quartets = line.split("\t")	### split quartet (4 number per position)
		for quartet in quartets:	### cycle through each quartet
			quad = (quartet.split())
			quadint = []
                        for basequal in quad:	### each intensity/number for each position
				quadint.append(int(basequal))
                        quadint.sort()
			quadint.reverse()
			basediff = quadint[0] - quadint[1]
                        #print "T=%i D=%i" % (quadint[0],basediff)

			if quadint[0] < threshold or basediff < difference:
				concat += "x"
			else:
				concat += "-"

  		head_match_regex = re.compile("\-{%i,%i}" % (consecutive,read_length)) 
		head_match = head_match_regex.search(concat)
 		if head_match != None:
			ok_read += 1
			col = head_match.span()
                        if not trim_info.has_key(read_number):
                                trim_info[read_number] = {}

			start = int(col[0])	
			end = int(col[1])
	
			trim_info[read_number]['start'] = start
			trim_info[read_number]['end'] = end

			if verbose:
				sub = concat[trim_info[read_number]['start']:trim_info[read_number]['end']]
				print "passed seqs:%i line#%i %s (start trim:%i,length:%i) %s\n" % (ok_read, read_number, concat, start, end, sub)

	LOG.write("%i out of %i sequences passed your filter (I >= %i and D >= %i and L >= %i)\n" % (ok_read, read_number, threshold, difference, consecutive))

	return trim_info


#--------------------------------------------------------------------------------------
def readNTrim(trim_info, seq, verbose, FASTA, LOG):

	"""         
        Parse a solexa/illumina sequence file and trim DNA sequence based user-defined intensity threshold/differences
	"""	
       

	read_number = 0
	gDNAlinker_count = 0
	usable_reads = 0

        dna_sequence_field = re.compile('^[ACTG]+$')
	gDNAlinker1_field = re.compile('^ATCCCC[GA]A')
	gDNAlinker2_field = re.compile('^ATCTAACAG')	

	if verbose:
		print "Printing trimmed sequences for all reads passing the set threshold values minus, excluding sequence containing linkers...\n"

        for line in seq:
		read_number += 1            ### tracks read number / will match order in prb file
        	line = line.rstrip('\r\n')
		info = line.split("\t")     ### split line, the seq file lists: lane tile xcoord y coord DNAseq 
		dna_string = info[4]
	
		if trim_info.has_key(read_number):
			trim_seq = dna_string[trim_info[read_number]['start']:trim_info[read_number]['end']]
			if re.match(dna_sequence_field, trim_seq):		### no ambiguous bases?
				if re.match(gDNAlinker1_field, trim_seq) or re.match(gDNAlinker2_field,trim_seq):	### matches gDNA linker?
					gDNAlinker_count += 1
				else:
					usable_reads += 1
					FASTA.write(">%s-%s-%s-%s\n%s\n" % (info[0],info[1],info[2],info[3],trim_seq))
					if verbose:
						print "line#%i %s (start trim:%i,length:%i) %s" % (read_number,info[4],trim_info[read_number]['start'],trim_info[read_number]['end'],trim_seq)
	LOG.write("%i out of %i sequences appear to be usable, after filtering out sequences hard-coded in this program * %i gDNA linker sequences*\n" % (usable_reads, read_number,gDNAlinker_count))
	return

if __name__ == '__main__':
	main()
	import time
	sys.exit()
