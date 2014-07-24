#!/usr/bin/env python

__doc__ = """
TQS

Trim Quality Sequences (TQS)

SYNOPSIS
   Quality trim FASTQ sequence reads using user-defined thresholds 
"""
__author__ = "Rene L. Warren"
__version__ = 'fastq'

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

#   Modified by Lance Parsons at Princton University's Lewis-Sigler Institute for Integrative Genomics
#   Adapted to trim "standard" FASTQ files (PHRED+33)

import sys, os, re, string, math
from datetime import datetime
from optparse import OptionParser


def main():
	usage = "Usage: %s --help"

	parser = OptionParser()
	parser.add_option("-f", "--fastq file", dest="fastqfile",
	                  help="fastq (fq) file - standard (ASCII+33) encoded PHRED quality scores / illumina (ASCII+64) encoded PHRED quality scores",)
        parser.add_option("-t", "--Phred quality threshold", dest="threshold", type="int", default=10,
                          help="Base intensity threshold value (Phred quality scores 0 to 40, default=10)",)
        parser.add_option("-c", "--consec", dest="consec", type="int", default=20,
                          help="Minimum number of consecutive bases passing threshold values (default=20)",)
        parser.add_option("-e", "--ASCII encoding type: 33 or 64", dest="encoding", type="int", default=64,
                          help="Type of ASCII encoding: 33 (standard) or 64 (illumina)  (default=64)",)
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
	                  help="Runs in Verbose mode.",)
	(opts, args) = parser.parse_args()
	
	try:
		f = open(opts.fastqfile)
		seq = f.readlines()
		f.close()
	except Exception, e:
		print "ERROR: Could not read from %s: %s" % (opts.fastqfile, e)
		print usage % (sys.argv[0:])
		sys.exit()


	fasta = "%s_T%sC%sE%s.trim.fa" % (opts.fastqfile, opts.threshold, opts.consec, opts.encoding)
	log = "%s.log" % opts.fastqfile
        minimum_length = 15


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
	
	if opts.consec < minimum_length:
		print "ERROR: -c must be a number larger than %i." % (minimum_length)
		sys.exit()

        if opts.encoding != 33 and opts.encoding != 64:
                print "ERROR: -e must be either 33 or 64."
                sys.exit()

	LOG.write("""
Running:
%s
-f %s
-c %s
-t %s
-e %s
Fasta file: %s

""" % (sys.argv[0:], opts.fastqfile, opts.consec, opts.threshold, opts.encoding, fasta))
	
        t1 = datetime.now()
        LOG.write("\n\nTrimming low quality bases: %s\n" % str(t1)[:len('2006-10-05 23:04')])
	readNtrim(seq, opts.threshold, opts.consec, opts.encoding, opts.verbose, FASTA, LOG)
        LOG.write("DNA sequences have been trimmed accordingly and placed in %s" % fasta)
	
	LOG.close()
	FASTA.close()
	return	

#--------------------------------------------------------------------------------------
def readNtrim(fastq, threshold, consecutive, encoding, verbose, FASTA, LOG):
	"""
	Return a Dictionary of sequence order number, with the index value and length to extract 
	"""
	trim_info = {}
	ok_read = 0
	read_number = 0
	record_line = 0
	
	if verbose:
		print "Printing trimming pattern for all reads passing the set threshold values...\n"
	
	for line in fastq:
		record_line += 1
		if record_line == 1:
			read_id = line.strip()
		elif record_line == 2:
			seq = line.strip()
		elif record_line == 3:
			qual_id = line.strip()
		elif record_line == 4:
			record_line = 0
			qual = line.strip()
			read_number += 1
			concat = ""			### concat builds a string of bases passing the user-defined filter 
			"""
			print "line%s\tseq:%s\tqual:%s\n" % (line,info[6],info[7])
			"""
			pos = 0
			for qual_char in qual:
				Q = (ord(qual_char) - encoding)
				pos += 1
				if Q < threshold:
					concat += "x"
				else:
					concat += "-"
				"""
				print "base#%i. Illumina qual (%s) == phredQ (%i)\n" % (pos,illumina_qual,Q)
				"""
	
			seq_len = len(seq)
	  		head_match_regex = re.compile("\-{%i,%i}" % (consecutive, seq_len)) 
			head_match = head_match_regex.search(concat)
	 		if head_match != None:
				ok_read += 1
				col = head_match.span()
	                        if not trim_info.has_key(read_number):
	                                trim_info[read_number] = {}
	
				start = int(col[0])	
				end = int(col[1])
	
	                        trim_seq = seq[start:end]
	                        FASTA.write(">%s\n%s\n" % (read_id, trim_seq))
	
				if verbose:
					print "%s\n%s\n%s\n passed seqs:%i line#%i %s (start trim:%i,end trim:%i) %s\n" % (read_id,seq,qual,ok_read, read_number, concat, start, end, trim_seq)

	LOG.write("%i out of %i sequences passed your filter (-t >= %i and -c >= %i)\n" % (ok_read, read_number, threshold, consecutive))

	return



if __name__ == '__main__':
	main()
	import time
	sys.exit()
