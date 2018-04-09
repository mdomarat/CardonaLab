#!/usr/bin/env python
import argparse
import ConfigParser

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from ssw_wrap import Aligner

#
# set on startup based on config file and command line arguments
#
# transposon sequence used to filter and strip reads
tn_seq = None
# reads must match tn_seq with a score better than this
min_score = None
# minimum length of match to tn_seq
min_match_length = None
# how much sequence after the tn is required
min_remaining_length = None


def set_params(args):
    global tn_seq
    global min_score
    global min_match_length
    global min_remaining_length
    global log
    global filtered_filename

    config = ConfigParser.SafeConfigParser()
    config.read(args.config)
    tn_seq=config.get("Defaults","tn_seq")
    min_score=config.getfloat("Defaults","min_score")
    min_remaining_length=config.getint("Defaults","min_remaining_length")
    min_match_length=config.getint("Defaults","min_match_length")
    filtered_filename=args.output

def filter_reads(readfile):
    print("Filtering reads\n")
    ssw = Aligner(tn_seq)
    total=0
    matched=0
    with open(filtered_filename,'w') as f:
        for title, seq, qual in FastqGeneralIterator(open(readfile)):
            total+=1
            res = ssw.align(seq,min_score, min_match_length)
            if res:
                end = res.query_end+1
                if len(seq)-end >= min_remaining_length:
                    matched+=1
                    f.write('@%s\n%s\n+\n%s\n' % (title, seq[end:], qual[end:]))
    print("%s of %s read had the tn seq\n" % (matched, total))


def main(argV=None):
    parser = argparse.ArgumentParser(description="filter out tn seq from reads")
    parser.add_argument('-c','--config', required=True,
        help='config file')

    parser.add_argument('-r','--reads', required=True,
        help='read file in fastq format')

    parser.add_argument('-o','--output', required=True,
        help='file to write filtered reads')

    args = parser.parse_args()
    set_params(args)
    filter_reads(args.reads)

if __name__ == "__main__":
    main()
