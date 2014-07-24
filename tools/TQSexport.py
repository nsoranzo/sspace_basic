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

import sys, os, re, string, math
from datetime import datetime
from optparse import OptionParser


def main():
	usage = "Usage: %s --help"

	parser = OptionParser()
	parser.add_option("-f", "--export file", dest="exportfile",
	                  help="Illumina export file - Output format from the Genome Analyzer",)
        parser.add_option("-t", "--Phred quality threshold", dest="threshold", type="int", default=10,
                          help="Base intensity threshold value (Phred quality scores 0 to 40, default=10)",)
        parser.add_option("-c", "--consec", dest="consec", type="int", default=20,
                          help="Minimum number of consecutive bases passing threshold values (default=20)",)
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
	                  help="Runs in Verbose mode.",)
	(opts, args) = parser.parse_args()
	
	try:
		f = open(opts.exportfile)
		seq = f.readlines()
		f.close()
	except Exception, e:
		print "ERROR: Could not read from %s: %s" % (opts.exportfile, e)
		print usage % (sys.argv[0:])
		sys.exit()


	fasta = "%s_T%sC%s.trim.fa" % (opts.exportfile,opts.threshold,opts.consec)
	log = "%s.log" % opts.exportfile
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

	LOG.write("""
Running:
%s
-f %s
-c %s
-t %s
Fasta file: %s

""" % (sys.argv[0:],opts.exportfile, opts.consec, opts.threshold, fasta))
	
        t1 = datetime.now()
        LOG.write("\n\nTrimming low quality bases: %s\n" % str(t1)[:len('2006-10-05 23:04')])
	readNtrim(seq, opts.threshold, opts.consec, opts.verbose, FASTA, LOG)
        LOG.write("DNA sequences have been trimmed accordingly and placed in %s" % fasta)
	
	LOG.close()
	FASTA.close()
	return	

#--------------------------------------------------------------------------------------
def readNtrim(export, threshold, consecutive, verbose, FASTA, LOG):
	"""
	Parse a solexa-illumina export file
	SOLEXA3_77_30V9CAAXX		4	1	1068	522		1	GGACAGCTGACAGCTGTTAAGAAGGACCCTATGTTAAAGGAAATGGATAC	YYYYYYYYYYYJYY
YYYYRYYYYYYYYYYYTTTTTOOOMOOOMMOOOOOG	chr13		36311743	F	50	52	121			187	R	N
	Return a Dictionary of sequence order number, with the index value and length to extract 
	"""
	trim_info = {}
	ok_read = 0
	read_number = 0

	if verbose:
		print "Printing trimming pattern for all reads passing the set threshold values...\n"
	
	for line in export:
		read_number += 1
		concat = ""			### concat builds a string of bases passing the user-defined filter 
		info = line.split() 	        ### split info 
		illumina_encoded_qual = list(info[7])
		"""
		print "line%s\tseq:%s\tqual:%s\n" % (line,info[6],info[7])
		"""
		pos = 0
		for illumina_qual in illumina_encoded_qual:
			pos += 1
			Q = 10 * math.log(1 + 10 ** ((ord(illumina_qual) - 64) / 10.0)) / math.log(10)
			if Q < threshold:
				concat += "x"
			else:
				concat += "-"
			"""
			print "base#%i. Illumina qual (%s) == phredQ (%i)\n" % (pos,illumina_qual,Q)
			"""

		seq_len = len(info[6])
  		head_match_regex = re.compile("\-{%i,%i}" % (consecutive, seq_len)) 
		head_match = head_match_regex.search(concat)
 		if head_match != None:
			ok_read += 1
			col = head_match.span()
                        if not trim_info.has_key(read_number):
                                trim_info[read_number] = {}

			start = int(col[0])	
			end = int(col[1])

			pair = ""
			if info[5] == "1":
				pair = "a"
			elif info[5] == "2":
				pair = "b"

                        trim_seq = info[6][start:end]
                        FASTA.write(">%s-%s-%s-%s%s\n%s\n" % (info[1],info[2],info[3],info[4],pair,trim_seq))

			if verbose:
				print "passed seqs:%i line#%i %s (start trim:%i,end trim:%i) %s\n" % (ok_read, read_number, concat, start, end, trim_seq)

	LOG.write("%i out of %i sequences passed your filter (-t >= %i and -c >= %i)\n" % (ok_read, read_number, threshold, consecutive))

	return



if __name__ == '__main__':
	main()
	import time
	sys.exit()
