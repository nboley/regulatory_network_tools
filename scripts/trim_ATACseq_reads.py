# Authors: Nathan Boley, Jason Buenrostro (Stanford University)
# Jason Buenrostro wrote the initial version, and Nathan Boley re-wrote
# to be more performant and use a more comprehensive mathcing system

##### IMPORT MODULES #####
# import necessary for python
import os
import re
import sys
import string
import itertools
import gzip

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import difflib

from argparse import ArgumentParser

import pyximport; pyximport.install()
from fuzz_align import fuzz_align

complement = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(sequence):
    return sequence.upper().translate(complement)[::-1]

def calc_edit_distance(l, r):
    return sum(x != y for x, y in itertools.izip(l, r))

# Align with mismatch, find first and move on, assumes only one
def fuzz_align_old(left_seq, right_seq, adapter, mismatch):
    """Find the optimal alignment of left_seq and right seq, 
    where the right portion of left seq is allowed to overlap
    the left portion of right seq. 

    """
    # loop through all allowable offsets
    rvs = []
    for offset in xrange(1, max(1, len(left_seq) - MIN_MATCH_LENGTH)): 
        l_primer = reverse_complement(left_seq[:offset])
        l_subset = left_seq[offset:offset+len(right_seq)]
        r_primer = right_seq[-offset:]
        dist = ( 
            calc_edit_distance(l_subset, right_seq)
            + calc_edit_distance(l_primer, adapter)/2.
            + calc_edit_distance(r_primer, adapter)/2. )
        if dist <= mismatch:  # find first then break
            rvs.append((offset, dist))
            
            
    # we didn't find a match
    if len(rvs) == 0: return None
    rvs.sort(key=lambda x:x[1])
    return rvs[0]

def parse_arguments():
    usage = "-a read1.fastq(.gz) -b read2.fastq(.gz) [--uncompressed] "
    opts = ArgumentParser(usage=usage)
    opts.add_argument("-a", required=True, 
                    help="<Read1> Accepts fastq or fastq.gz")
    opts.add_argument("-b", 
                    help="<Read2> Accepts fastq or fastq.gz")
    opts.add_argument("--uncompressed", "-u", action="store_true", 
                    default=False, help="Print uncompressed output file")
    options = opts.parse_args()
    
    # name input and outputs
    p1_in = options.a
    if options.b == None:
        p2_in = options.a.replace("R1", "R2")
    else:
        p2_in = options.b
    
    # name outputs and print to working dir
    p1_file = p1_in.split('/')[-1]
    p2_file = p2_in.split('/')[-1]

    #check for file type and open input file
    append = p1_in.split('.')[-1]
    if append == "fastq":
        p1_rds = open(p1_in,'r')
        p2_rds = open(p2_in,'r')
        p1_out = re.sub(".fastq", ".trim.fastq", p1_file)
        p2_out = re.sub(".fastq", ".trim.fastq", p2_file)
    elif append == "fq":
        p1_rds = open(p1_in,'r')
        p2_rds = open(p2_in,'r')
        p1_out = re.sub(".fq", ".trim.fastq", p1_file)
        p2_out = re.sub(".fq", ".trim.fastq", p2_file)
    elif append == "gz":
        p1_rds = gzip.open(p1_in,'r')
        p2_rds = gzip.open(p2_in,'r')
        p1_out = re.sub(".fastq.gz", ".trim.fastq", p1_file)
        p2_out = re.sub(".fastq.gz", ".trim.fastq", p2_file)
    else:
        sys.exit("ERROR! The input file2 must be a .fastq or .fastq.gz")

    # initialize write files
    if options.uncompressed == False:
        print "WARNING: GZIP output can be very slow (it's the python gzip module's fault)"
        r1_write = gzip.open(p1_out+'.gz', 'wb')
        r2_write = gzip.open(p2_out+'.gz', 'wb')
    elif options.uncompressed == True:
        r1_write = open(p1_out, 'w')
        r2_write = open(p2_out, 'w')
    
    return ( 
        (FastqGeneralIterator(p1_rds), FastqGeneralIterator(p2_rds)), 
        (r1_write, r2_write) )

def get_next_read(p1_rds, p2_rds):
    for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in itertools.izip(
            p1_rds, p2_rds):
        assert f_id.split()[0] == r_id.split()[0], (
            "Read ids do not match (%s, %s)" % (f_id, r_id))
        yield (f_id, r_id), (f_seq, r_seq), (f_q, r_q)

def main():
    print >> sys.stderr, "="*80
    print >> sys.stderr, "THIS IS NOT PRODUCTION CODE AND HAS BUGS - USE AT YOUR OWN RISK"
    print >> sys.stderr, "THE FOLLOWING OPTIONS ARE HARD CODED:"
    # the maximum number of mismatches we allow for aligned sequences
    mismatch=10  
    print >> sys.stderr, "MAX NUM MISMATCH: ", mismatch
    adapter = "CTGTCTCTTATACACATCT" # illumina nextera adapter
    print >> sys.stderr, "ADAPTER SEQUENCE (illumina nextera adapter): ", adapter
    # the minimum required pair overlap for reads to be considered overlapping
    min_match_length = 24 
    print >> sys.stderr, "MIN SEQUENCE MATCH LENGTH: ", min_match_length
    print >> sys.stderr, "="*80

    (p1_rds, p2_rds), (p1_ofp, p2_ofp) = parse_arguments()


    """
    # debugging code
    seq1 = "AAA" + adapter
    seq2 = "TTT" + adapter
    align_1 = fuzz_align(
        seq1, seq2, adapter, mismatch, MIN_MATCH_LENGTH) 
    align_2 = fuzz_align_old(
        seq1, seq2, adapter, mismatch) 
    print align_1, align_2
    print seq1
    print seq2

    try: idx, mis = align_2
    except: idx, mis = align_1

    print adapter
    print seq2[-idx:]

    print adapter
    print seq1[-idx:]

    print reverse_complement(seq1[:-idx])
    print seq2[:-idx]
    """
    #return
    for i, ((seqhead1, seqhead2), (seq1, seq2), (qual1, qual2)) in enumerate(
                get_next_read(p1_rds, p2_rds)):
        #if i > 500: break
        if i%100000 == 0: print i
        
        align_1 = fuzz_align(
            seq1, seq2, adapter, mismatch, MIN_MATCH_LENGTH) 
        ## Debugging code
        #align_2 = fuzz_align_old(
        #    rc_seq1, seq2, adapter, mismatch) 
        
        if align_1 != None:
            #print seq1
            #print seq2
            idx, mis = align_1
            #except: idx, mis = align_1
            #print align_1, align_2
            
            #print adapter
            #print seq2[-idx:]

            #print adapter
            #print seq1[-idx:]
            if idx > 0:
                seq1 = seq1[:-idx]
                qual1 = qual1[:-idx]

                seq2 = seq2[:-idx]
                qual2 = qual2[:-idx]
            #print reverse_complement(seq1)
            #print seq2
        
        p1_ofp.write("@%s\n%s\n+\n%s\n" % (
            seqhead1, seq1, qual1 ) )

        p2_ofp.write("@%s\n%s\n+\n%s\n" % (
            seqhead2, seq2, qual2 ) )

if __name__ == '__main__':
    main()
