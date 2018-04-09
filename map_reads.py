#!/usr/bin/env python
# determine transposon insertion sites.
import argparse
import ConfigParser
import csv
import os

import re
import subprocess
from itertools import groupby
import random

from Bio import SeqIO

def bowtie2_build(ref,outbase):
    output = subprocess.check_output(["bowtie2-build", ref, outbase])
    return output

def bowtie2(readfile, reference, outfile):
     try:
       output = subprocess.check_output(["bowtie2", "--end-to-end", "-p 4", "-a", "-x", reference, "-U", readfile, "-S", outfile], stderr=subprocess.STDOUT)
     except subprocess.CalledProcessError as e:
       print e.output
     return output

def read_samfile(infile):
    with open(infile,'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0][0] == '@':
                # skip header rows
                continue
            else:
                yield row

#should nicely separate CIGAR entries
cigar_pat = re.compile(r"\d+[MIDNSHP=X]{1}")
# cigar2end adapted from: https://www.biostars.org/p/48333/

def cigar2end( left,cigar ):
  """Return right-most position of aligned read."""
  #store info about each CIGAR category
  counts={ "M":0, #M 0 alignment match (can be a sequence match or mismatch)
           "I":0, #I 1 insertion to the reference
           "D":0, #D 2 deletion from the reference
           "N":0, #N 3 skipped region from the reference
           "S":0, #S 4 soft clipping (clipped sequences present in SEQ)
           "H":0, #H 5 hard clipping (clipped sequences NOT present in SEQ)
           "P":0, #P 6 padding (silent deletion from padded reference)
           "=":0, #= 7 sequence match
           "X":0, #X 8 sequence mismatch
        }
  #split cigar entries
  for centry in cigar_pat.findall(cigar):
    ccount  = int(centry[:-1])
    csymbol = centry[-1]
    counts[csymbol] += ccount
  #get number of aligned 'reference' bases
  aligned = counts["M"] + counts["D"] + counts["N"] + counts["="] + counts["X"]
  right   = left + aligned
  return right


def get_percentage ( init_string, backwards):
  MAX_LENG = 25

  dir = 1
  caret_posn = 0 
  string = init_string
  regex = r'^([A-Z]|\^[A-Z]+)(\d+)(.*)'
  if backwards:
    dir = -1 
    caret_posn = -1
    string = init_string[::dir]
    regex = r'^([A-Z]|[A-Z]+\^)(\d+)(.*)'
 
  matc = re.match(r'^(\d+)(.*)$', string)
  total_size = int(matc.group(1)[::dir])
  total_mismatch = 0
  rest = matc.group(2)
  
  while rest and total_size < MAX_LENG:
     matc = re.match(regex, rest)
     mismatch = matc.group(1)  
     if mismatch[caret_posn] == '^':
       mismatch_leng = len(mismatch)-1
     else:
       mismatch_leng = 1
     total_mismatch += mismatch_leng
     total_size += mismatch_leng

     # we've gone over with mismatches, subtract them.
     if total_size > MAX_LENG:
       diff = total_size - MAX_LENG
       total_mismatch -= diff 
       total_size -= diff
       rest = None
     else: # we need more, try matches     
       match_leng = int (matc.group(2)[::dir] )
       if len(matc.groups()) > 2:
         rest = matc.group(3)
       else:
         rest = None
       total_size += match_leng
       
       if total_size > MAX_LENG:
         total_size = MAX_LENG  # gone over with matches.

  percentage = float( total_size - total_mismatch ) / total_size  
  return percentage 

# get the first/last 25 positions (pctg)   
def parse_md_z_25 (md_z_string):
  front_percentage = get_percentage ( md_z_string, False) 
  back_percentage = get_percentage ( md_z_string, True)
  return (front_percentage, back_percentage)
  

# get the edges of the MD:Z string -- all the matches
# EXCEPT the 0s -- they are ignored.
def parse_md_z (md_z_string):
  # look at edges
  front_zero = re.match(r'^0[A-Z](.*)$', md_z_string)
  md_z_string_F = md_z_string
  while front_zero:
     md_z_string_F = front_zero.group(1) 
     front_zero = re.match(r'^0[A-Z](.*)$', md_z_string_F)
  
  front_match = re.match(r'^(\d+)' , md_z_string_F)
  front_size = -1
  if front_match:
     front_size  = int(front_match.group(1))

  back_zero = re.match(r'^(.*)[A-Z]0$', md_z_string)
  md_z_string_B = md_z_string
  while back_zero:
     md_z_string_B = back_zero.group(1) 
     back_zero = re.match(r'^(.*)[A-Z]0$', md_z_string_B)

  back_match = re.search(r'(\d+)$' , md_z_string_B)
  back_size = -1
  if back_match:
     back_size  = int(back_match.group(1))
  

  return [front_size,back_size]


# filters rows by % of first 25 bp matching. 
def filter_by_percentage( rows, pctg):
  for row in rows:
       flag = int(row[1])
       md_match = None
       ind = -1
       while (not md_match) and ind >= -1*len(row):
         md_match = re.match(r'MD:Z:(\S*)', row[ind])
         ind -= 1
       md_string = 'empty'
       if md_match:
         md_string = md_match.group(1)
         pctgs = parse_md_z_25 (md_string)

       if (pctgs[flag&16==16]*100. >= pctg): 
         yield row

      
# filters rows by length of initial match  
def filter_short_matches(rows, leng):
    for row in rows:
       flag = int(row[1])
       md_match = None
       ind = -1
       while (not md_match) and ind >= -1*len(row):
         md_match = re.match(r'MD:Z:(\S*)', row[ind])
         ind -= 1
       md_string = 'empty'
       if md_match:
         md_string = md_match.group(1)
         match_lengths = parse_md_z (md_string)
      
       if ((flag&16==16 and match_lengths[1] >= leng) or
           (flag&16==0  and match_lengths[0] >= leng)):
         yield row

# checks sam output to see if flag&4 is unset, i.e., no matches
def filter_unmaped(rows):
    for row in rows:
        if not int(row[1])& 4:
            yield row


def filter_duplicates(rows):
    for k, g in groupby(rows,lambda x: x[0]):
        # need to check for differences in alignment quality
        g = list(g)
        yield random.choice(g)


def alignments_to_insertions(rows):
    counts = {} # dict ()

    for row in rows:
        src = row[0]
        ref = row[2]  # source -- chromosome, etc. 

        ref_counts = counts.get(ref,{}) # dict() is default
        if not ref_counts:
            counts[ref] = ref_counts

        lpos = int(row[3]) # leftmost position on forward strand 
        rpos = lpos #default if forward strand
        flag = int(row[1])
        if (flag & 16 ) == 16:
          rpos = cigar2end(rpos,row[5]) 
        c = ref_counts.get(rpos,[0,0])
        if c == [0,0]:
            ref_counts[rpos] = c

        if (flag & 16) == 16:
            c[1] += 1
        else:
            c[0] += 1

    return counts

def write_insertions(fname,counts):
    with open(fname,'wb') as f:
        f.write("reference\tposition\tstrand\tcount\n")
        for ref,ref_counts in sorted(counts.items()):
            for rpos,v in sorted(ref_counts.items()):
                if v[0] > 0:
                    f.write("%s\t%d\t%d\t%d\n" % (ref,rpos,1,v[0]))
                if v[1] > 0:
                    f.write("%s\t%d\t%d\t%d\n" % (ref,rpos,-1,v[1]))

def map_insertions(readfilename, genomefilename, outfilename, length, pctg):
    fastafilename="test.fasta"
    samfilename= os.path.basename(outfilename) + '.sam'   # was "mapped.sam"
    SeqIO.convert(genomefilename,"genbank",fastafilename,"fasta")
    print bowtie2_build(fastafilename,"ref")
    print "Mapping reads:"
    print bowtie2(readfilename,"ref",samfilename)
    print "Filtering reads"
    counts = alignments_to_insertions(
        filter_by_percentage(
	filter_short_matches(
        filter_duplicates(
            filter_unmaped(
                read_samfile(samfilename))),length), pctg ))
    return counts





def main(argV=None):
    parser = argparse.ArgumentParser(description="determine tn insertion sites")

    parser.add_argument('-r', '--reads', required=True,
        help='read file in fastq format')

    parser.add_argument('-g', '--genome', required=True,
        help='reference genome in genbank format')

    parser.add_argument('-o', '--output', required=True,
        help='file to write insertions')

    parser.add_argument('-l', '--length', type=int, required=False, default=0, 
        help='minimum length of initial match. if absent, all matches reported')

    parser.add_argument('-p', '--percentage', type=int, required=False, default=0, 
        help='minimum percentage agreement in first 25 bp. of match. if absent, 0.')

    args = parser.parse_args()

    insertions = map_insertions(args.reads, args.genome, args.output, args.length, args.percentage)
    write_insertions(args.output, insertions)

if __name__ == "__main__":
    main()
