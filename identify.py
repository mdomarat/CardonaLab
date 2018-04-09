#!/usr/bin/env python
# identify genes that are downstream from insertion sites,
# based on certain criteria provided by the user
# mdomarat, june 2016.
# 
# updated july 2016 to incorporate margin computations.
# update, july 20, 2016 - fix bugs that duplicate lines. 
# update, aug 23, 2016 -- removed strand mismatch error. raised due to reannotation. 
# update, jan 22, 2018 -- fixed bug: when looking backwards in elif, didn't use corrected start,end
import argparse
import ConfigParser
import csv
import sys

import re
import subprocess
from collections import defaultdict
from itertools import groupby, tee, izip_longest,chain,izip
import random

from Bio import SeqIO

qualifiers_of_interest = ['locus_tag','old_locus_tag','protein_id']
def check_op ( op_list ,  gene_info):
   retVal = []
   posn = int ( gene_info [3] )
   orient = -1 if ( gene_info [ 5] == '-' ) else 1
   start_ind = 3 if orient==1 else 4
   #print 'orient',orient, start_ind
   for x in op_list:
      if ( orient*int(gene_info [ start_ind ]) <  orient*int(x [ start_ind ]) ):
        #print '-->', x, gene_info [ start_ind]*orient, orient*x[start_ind]
	retVal.extend( x[2:] )
   return retVal
   
def add_operon_info(rows,operon_filename):
    if operon_filename:
	    dicts = get_op_dict(operon_filename)
	    operons_by_bcal = dicts[0] # to find if the current cds has a bcal that's in the operon file.
	    operons_by_id  = dicts[1] 
	    for row in rows:
	       if (len(row) >= 9 and row[8] in operons_by_bcal): ## is there a BCAL AND the BCAL is in the operon file?
		 op_info = operons_by_bcal[row[8]] 
		 row.extend ( op_info [6:9] ) ## add length, cog number, product for FIRST gene
		 #if ((row[2]=='1' and op_info[5] == '+') or (row[2] == '-1' and op_info[5] == '-')):
		 #  raise Exception( 'strand mismatch error')
		 if ( len( operons_by_id [ op_info[0] ] ) > 1 ) :
	 
		   row.extend(check_op (operons_by_id[ op_info[0] ], op_info ))
	       yield row
    else:
	for row in rows: yield row


def get_op_dict(infile):
  ops_by_bcal = dict() 
  ops_by_id = defaultdict(list)
  with open(infile,'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if (len(row) == 0 or row[0][0] == 'O'):
                # skip header rows that start with "O"peron
                continue
            else:
              ops_by_bcal [ row[2] ] = row;
              ops_by_id   [ row[0] ].append(row)
  return (ops_by_bcal, ops_by_id ) 


## this looks for the next cds in the direction indicated by dir 
## and sees if it is within the maxDist. (if maxDist is really big, this will run a long time 
## if it doesn't find anything).
## note: should only be used when we are trying to find the next CDS, not the current one, obv.
## dir should be +1/-1
def next (cds_list, dir, ind, maxDist, posnAtIns):
  curr = ind + dir ## find the next one
  found = False
  tooFar = False
  marker = -1 if dir < 0 else len(cds_list) 
  while (curr != marker and not found and not tooFar):
    curr_cds = cds_list[curr]
    curr_bound = [ curr_cds.location.start, curr_cds.location.end]
    curr_close_end = curr_bound [ (1 - dir)/2 ] 
    curr_far_end = curr_bound [ (dir + 1)/2 ]  ## needed?? 
    if (dir == curr_cds.strand and dir*posnAtIns <  dir*curr_close_end and dir*(posnAtIns+maxDist) < dir*curr_close_end):
       found = True
    elif ( dir == curr_cds.strand and dir*(posnAtIns+maxDist) >= dir*curr_close_end):
       tooFar = True
    else:
       curr += dir 
  if found:
    return cds_list[curr]
  else:
    return None

def read_tsv(infile,genomefilename,leng,left,right):
    cds_dict = get_CDS_dict(genomefilename)
    with open(infile,'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
          if row[0][0] == 'r':
                # skip header rows
                continue
	  else:
               all_cds = cds_dict[row[0]] # get CDS on specified seq 
               posn_ins = int ( row[1] )  # where is the insertion position?
               strand_ins = int (row [2] )  # which strand is the insertion on?
	       found = False
               know_bad = False
               distance_reversed = sys.maxint # have we found the closest reverse strand CDS?
               r_extend = None  ## information to put on the end of the row if we have a closest rev strand CDS
	       for ind in range(len(all_cds)):
                 cds = all_cds[ind] ## current cds under consideration.
                 if strand_ins != cds.strand: ## insertion should be on the opposite strand 
		   bound = (cds.location.start, cds.location.end)
                   if (cds.strand == -1):
 			bound = bound[::-1]
                   start,end = bound

                   if (cds.strand*start <= cds.strand*posn_ins and cds.strand*end >= cds.strand*posn_ins):
                     ## insertion is INSIDE the CDS (with correct orientation) -- check margins
                     left_margin_distance = left*len(cds)/100
                     right_margin_distance = right*len(cds)/100
                     if (not know_bad and cds.strand*(posn_ins - start )  < left_margin_distance):
                       row.extend([str(cds.strand),str(cds.location.start+1), str(cds.location.end)])
                       for x in qualifiers_of_interest:
                        if x in cds.qualifiers:
				row.append(cds.qualifiers[x][0]  )
			else:
				row.append('N/A')
                       if 'db_xref' in cds.qualifiers:
                     	row.extend(cds.qualifiers['db_xref'])
                       found = True
                       break
                     elif (not know_bad and cds.strand*(end - posn_ins) < right_margin_distance): 
                       possible_cds = next(all_cds,cds.strand,ind,leng,posn_ins)
                       if possible_cds:
                         #print strand_ins,'end margin',posn_ins,possible_cds
                         row.extend([str(possible_cds.strand),str(possible_cds.location.start+1), str(possible_cds.location.end)])
                         for x in qualifiers_of_interest:
                           if x in possible_cds.qualifiers:
  				row.append(possible_cds.qualifiers[x][0]  )
      		           else:
  				row.append('N/A')
                         if 'db_xref' in possible_cds.qualifiers:
                     	   row.extend(possible_cds.qualifiers['db_xref'])
                         found = True
 			 break 
                     else:
                       # we've found something that's in the middle (not on the margins).  discard and continue.
                       row.append(str(cds.strand) + ' N/A - M') 
		       row.extend([str(cds.location.start+1), str(cds.location.end)])
                       for x in qualifiers_of_interest:
                         if x in cds.qualifiers:
                              row.append(cds.qualifiers[x][0]  )
                         else:
                              row.append('N/A')
                       if 'db_xref' in cds.qualifiers:
                         row.extend(cds.qualifiers['db_xref'])
                       
                       know_bad = True
                       found = True
                       break
                        
                   elif (not know_bad and cds.strand*start > cds.strand*posn_ins and cds.strand*(start - posn_ins  ) < leng):
                     ## insertion is close enough to the next one to add details 
                     ## note: bug fixed jan/18 -- wasn't using start/end still using old cds.location.*
                     if (cds.strand > 0):  ## forward strand -- easy. you've found the next one. 
                       row.extend([str(cds.strand),str(cds.location.start+1), str(cds.location.end)])
                       for x in qualifiers_of_interest:
                          if x in cds.qualifiers:
  				row.append(cds.qualifiers[x][0]  )
			  else:
				row.append('N/A')
                       if 'db_xref' in cds.qualifiers:
                     	 row.extend(cds.qualifiers['db_xref'])
                       found = True
                       break
                     else:
                       ## reverse strand. 
                       found = True
                       dist = abs ( cds.location.start - posn_ins ) ## double check in future - should use start?
                       if (dist < distance_reversed ):
                          r_extend = [ str(cds.strand), str(cds.location.start+1), str( cds.location.end) ] 
                          for x in qualifiers_of_interest:
                            if x in cds.qualifiers:
  				r_extend.append(cds.qualifiers[x][0]  )
			    else:
				r_extend.append('N/A')
                          if 'db_xref' in cds.qualifiers:
                     	     r_extend.extend(cds.qualifiers['db_xref'])
	       if not found and not know_bad:  ## nothing happened -- too far away, not in middle, etc.
                 row.append('N/A')
                 know_bad = True
               elif r_extend and not know_bad: ## we found a reverse strand, and now know the smallest.  
                 row.extend(r_extend)
               yield row
 

# get CDS dictionary -- entries should be sorted by "start" positions
def get_CDS_dict(genomefilename):
    d = defaultdict(list) 
    for seq_record in SeqIO.parse(genomefilename, "genbank"):
      for feat in seq_record.features:
        if feat.type == 'CDS':
          d [ seq_record.id ].append(feat);
    for x in d:
        d[x].sort(key=lambda cds:cds.location.start) # sort by starting positions 
    return d
 
def write_additional_info( rows, outfile):
    with open(outfile,'wb') as f:
      f.write("reference\tposition\tstrand\tcount\tCDS strand\tCDS start\tCDS end\tLocus TAG\tOld Locus Tag\tGI/Gene IDs\tlength\tCOG\tproduct\tOperon Info\n")
      for row in rows:
        f.write ( '\t'.join(row) )
        f.write ( '\n'  )

def main(argV=None):
    parser = argparse.ArgumentParser(description="identify downstream genes from insertion sites")

    parser.add_argument('-t', '--tsv', required=True,
        help='tsv file of insertion sites [output of map_reads.py]')

    parser.add_argument('-g', '--genome', required=True,
        help='reference genome in genbank format')

    parser.add_argument('-p', '--operon', required=False, 
        help='operon file in door tabular format')

    parser.add_argument('-o', '--output', required=True,
        help='file to write new tsv file')

    parser.add_argument('-l', '--length', default=sys.maxint, type=int, required=False,
        help='maximum downstream distance to gene (if absent, default none)')

    parser.add_argument('-L', '--leftmargin', default=20, type=int, required=False,
        help='maximum 5\' "margin" percentage (if absent, default 20%)')

    parser.add_argument('-R', '--rightmargin', default=20, type=int, required=False,
        help='maximum 3\' "margin" percentage (if absent, default 20%)')

    args = parser.parse_args()

    print ('starting...')
    x = add_operon_info(read_tsv(args.tsv,args.genome,args.length,args.leftmargin,args.rightmargin),args.operon)
    print ('done gathering data, writing.')
  
    write_additional_info( x , args.output )
    print ('done')

if __name__ == "__main__":
    main()

