#!/usr/bin/env python
# merge control and treatment CSV files 
# mdomarat, june 2016.
# january 2018 -- added nargs  

from __future__ import division
import argparse
import ConfigParser
import csv

import re
import subprocess
from collections import defaultdict
from itertools import groupby
import random

from Bio import SeqIO
import math


# merges two dictionaries d1 and d2
# dictionaries have values that are lists.
# d1's lists should have length leng1
# d2's lists should have length leng2
# the joined lists have length leng1+leng2
def merge_dicts (d1,d2,leng1,leng2):
  merged = dict()
  for x in d1:
    if x in d2:
      merged [ x] =  d1[x] + d2[x]  
    else:
      merged [ x] =  d1[x] + [0]*leng2 
  for x in d2:
     if not ( x in merged) :
         merged [ x ] =  [0]*leng1 + d2[x]   
  for x in merged:
    assert(len(merged[x]) == leng1+leng2)  
  return merged

# simpler marge
def merge_simple_dicts (d1,d2):
  merged = d1.copy() 
  merged.update(d2)
  return merged
     

# get insertion sites from a list of either
# control or treatment files
def get_ins_sites( f_list):
  ins_sites = dict()
  additional_info = dict()
  for counter,filename in enumerate(f_list):
    with open(filename,'rb') as f:
      reader = csv.reader(f, delimiter='\t')
      for row in reader:
        if (row[0] == 'reference'):
          continue
        else:
          index = (row[0], int(row[1]), int(row[2]))
          if (index in ins_sites):
             ins_sites[index].append(row[3])
          else:
             ins_sites[index] = [0]*counter + [int(row[3])]
             additional_info [ index ]  = row[4:] 
          #ins_sites [ ( row[0], int(row[1]), int(row[2]) ) ]  = row[3], filename
  for x in ins_sites:
    ins_sites[x].extend( [0]*( len(f_list) - len(ins_sites[x]) )) 
    assert( len(ins_sites[x]) == len(f_list) )
  return (ins_sites, additional_info)

def read_tsv(control_list, treatment_list,outfile):
    control_ins_sites, control_addl_info = get_ins_sites ( control_list) 
    treat_ins_sites, treat_addl_info = get_ins_sites( treatment_list)
    merged_ins_sites = merge_dicts ( control_ins_sites, treat_ins_sites, len(control_list), len(treatment_list))
    merged_addl_info = merge_simple_dicts ( control_addl_info, treat_addl_info) 
    file_labels = '\t'.join(control_list + treatment_list) 
    with open(outfile,'wb') as f:
      f.write("reference\tposition\tstrand\t" + file_labels + "\tCDS strand\tCDS start\tCDS end\tLocus TAG\tOld Locus Tag\tGI/Gene IDs\tDOOR info\n")
      for x in sorted(list(merged_ins_sites)):
        values = merged_ins_sites[x]
        out_string = '\t'.join ( ( x[0] , str(x[1]) ,  str(x[2])  )  ) + '\t' + '\t'.join([str(z) for z in merged_ins_sites[x]])
        f.write (out_string + '\t' ) 
        f.write ( '\t'.join ( merged_addl_info [ x] ) if len(merged_addl_info[x]) > 1 else merged_addl_info[x][0] ) 
        f.write ( '\n'  )

def main(argV=None):
    print('starting...')
    parser = argparse.ArgumentParser(description="merge")

    parser.add_argument('-c', '--control', required=True, nargs ='*',
        help='tsv file of control insertion sites')

    parser.add_argument('-t', '--treatment', required=True, nargs='*',
        help='tsv file of treatment insertion sites')

    parser.add_argument('-o', '--output', required=True,
        help='file to write new tsv file')

    args = parser.parse_args()

    read_tsv(args.control,args.treatment, args.output)
    print('done. output files is ' + args.output)
  

if __name__ == "__main__":
    main()
