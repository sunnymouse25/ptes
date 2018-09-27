# Takes segemehl.sam.nohead as input, makes junction table and HTML table for junctions counts

# Imports
import subprocess, sys, os, pickle, json, yaml
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np
from interval import interval

from ptes.lib.general import init_file, writeln_to_file, shell_call, write_to_file
from ptes.ptes import get_read_interval, one_interval, get_interval_length
from ptes.ucsc.ucsc import order_interval_list, list_to_dict, get_track_list

### Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="Segemehl output file; SAM file without header")
parser.add_argument("-j","--json", type=str,
                    default = 'no',
                    help="if json file with junctions from GTF is present, write yes")   
parser.add_argument("-a","--acceptors", type=str,
                    default = 'gtf_acceptors.json',
                    help="json file with GTF acceptors, optional")
parser.add_argument("-d","--donors", type=str,
                    default = 'gtf_donors.json',
                    help="json file with GTF donors, optional")                    
parser.add_argument("-o","--output", type=str,
                    help="Output folder for .csv and .html files")                                       
args = parser.parse_args()

### Functions

def read_coord(intervals, infos, read_name, cigar, leftpos, XI, XQ, chrom, flag):           
        xq = int(XQ)
        if flag & 16 == 0:   
            chain = '+'
        else:   
            chain = '-'
        infos[read_name][XI].append(chrom)
        infos[read_name][XI].append(chain)
        read_interval = one_interval(get_read_interval(cigar,leftpos, output = 'interval'))   # local alignment is always bound by M
        intervals[read_name][XI].append((xq,read_interval))
        
def BinarySearch(a, key, frm, to):
    if to-frm == 0: 
        return False
    i = (to+frm)//2 
    if a[i] == key:
        return True
    elif a[i] > key:
        return BinarySearch(a,key,frm,i)
    else:
        return BinarySearch(a,key,i+1,to) 

        
### Main  
segemehl_outfile = args.input  # SAM file without header
dirname = os.path.dirname(os.path.realpath(segemehl_outfile))
if dirname == '':
    dirname = '.'
path_to_file = args.output.rstrip('/')

col_names = ['read_name','flag','chrom','leftpos','cigar','xi','xq']
# col_nums = [0,1,2,3,5,14,19]
tag_list = ['XI','XQ']

read_intervals = defaultdict(lambda: defaultdict(list))   # mapped intervals
read_infos = defaultdict(lambda: defaultdict(list))   # mapped chrom(s) and chain(s)

# Reading SAM input

with open(segemehl_outfile, 'r') as df_segemehl:
    for line in df_segemehl:
        row = line.strip().split('\t')   
        read_name = row[0]
        flag = int(row[1])
        chrom = row[2]
        leftpos = row[3]
        cigar = row[5]
        sam_attrs = {'read_name' : read_name,
                        'flag': flag,
                        'chrom' : chrom,
                        'leftpos' : leftpos,
                        'cigar' : cigar}
        tags = dict.fromkeys(tag_list, None)
        for elm in row[14:]:
            for tag in tag_list:
                if tag in elm:
                    tags[tag] = elm
        sam_attrs.update(tags)       
        sam_attrs['XI'] = sam_attrs['XI'].strip('XI:i:')   # Mates have different xi!        
        if sam_attrs['XQ']:
            sam_attrs['XQ'] = sam_attrs['XQ'].strip('XQ:i:')
            read_coord(read_intervals, read_infos, **sam_attrs)
            
with open('%s/intervals_df_segemehl.json' % path_to_file, "w") as intervals_json,\
    open('%s/infos_df_segemehl.json' % path_to_file, "w") as infos_json:
    json.dump(read_intervals, intervals_json)
    json.dump(read_infos, infos_json)   
    
# Exons GTF to junctions dict
if args.json == 'no':
    gtf_exons_name = '/uge_mnt/home/sunnymouse/Human_ref/hg19_exons.gtf'
    gtf_donors = defaultdict(set)
    gtf_acceptors = defaultdict(set)
    with open(gtf_exons_name, 'r') as gtf_exons_file:
        for line in gtf_exons_file:
            line_list = line.strip().split()
            chrom = line_list[0]
            strt = int(line_list[3])
            end = int(line_list[4])
            chain = line_list[6]
            if chain == '+':
                gtf_donors[chrom].add(end+1)
                gtf_acceptors[chrom].add(strt-1)
            elif chain == '-':
                gtf_donors[chrom].add(strt-1)
                gtf_acceptors[chrom].add(end+1)
    
    gtf_donors_list = {}   # sets are not JSON serializable
    gtf_acceptors_list = {}
    
    for key in gtf_donors:
        gtf_donors_list[key] = list(gtf_donors[key])
    for key in gtf_acceptors:
        gtf_acceptors_list[key] = list(gtf_acceptors[key])
    
    with open('%s/%s' % (path_to_file, 'gtf_donors.json'), 'w') as donors_file,\
         open('%s/%s' % (path_to_file, 'gtf_acceptors.json'), 'w') as acceptors_file:
        json.dump(gtf_donors_list, donors_file)
        json.dump(gtf_acceptors_list, acceptors_file)
else:
    with open(args.donors, 'r') as donors_file,\
        open(args.acceptors, 'r') as acceptors_file:
        gtf_donors = yaml.safe_load(donors_file)
        gtf_acceptors = yaml.safe_load(acceptors_file)
    for key in gtf_donors:        
        gtf_donors[key] = map(int, gtf_donors[key])
        gtf_donors[key] = set(gtf_donors[key])
    for key in gtf_acceptors:
        gtf_acceptors[key] = map(int, gtf_acceptors[key])
        gtf_acceptors[key] = set(gtf_acceptors[key])
    
    
# Mapped junctions table            
intervals_list = []   # to remember read intervals
junc_list = []   # from mapped read intervals to list of junctions

for key in read_intervals.keys():   # key is read_name    
    donor_ss = np.nan
    acceptor_ss = np.nan
    xi = 0
    while read_intervals[key].get(str(xi), None):   # xi - number of current read alignment                
        tuples = read_intervals[key][str(xi)]   # list of intervals of current read alignment    
        tuples = sorted(tuples, key=lambda x:x[0])   # sort by xq        
        values =  [i[1] for i in tuples]   # get rid of xq
        values = order_interval_list(values)   # ascending order is essential for BED lines
        n_j = len(values) - 1
        infos = read_infos[key][str(xi)]
        chroms = set(infos[::2])
        chains = set(infos[1::2])
        xi += 1
        if len(chroms) == 1 and len(chains) == 1:    # read must be mapped to the same chrom and chain
            chrom = chroms.pop()
            if chrom.startswith('chr'):    # skip contigs
                chain = chains.pop()                    
                if n_j > 0:                 # should be every read, just to be sure                                        
                    for i, current_interval in enumerate(values):
                        if i > 0:
                            if chain == '+':                        
                                donor_ss = int(values[i-1][0].sup) + 1
                                acceptor_ss = int(values[i][0].inf) - 1
                            elif chain == '-':                        
                                donor_ss = int(values[i-1][0].inf) - 1
                                acceptor_ss = int(values[i][0].sup) + 1    
                            donors = gtf_donors[chrom]
                            acceptors = gtf_acceptors[chrom]                            
                            annot_donor = donor_ss in donors 
                            annot_acceptor = acceptor_ss in acceptors
                            chimeric = True if donor_ss > acceptor_ss else False
                            junc_list.append({'read_name' : key,
                                                               'aln' : xi-1,
                                                               'n_junctions' : n_j,                                
                                                               'chrom' : chrom, 
                                                               'chain' : chain,
                                                               'donor' : str(donor_ss),
                                                               'annot_donor' : annot_donor,
                                                               'acceptor' : str(acceptor_ss),
                                                               'annot_acceptor' : annot_acceptor,
                                                               'chimeric' : chimeric}) 

  
mapped_junc_df = pd.DataFrame(junc_list)
mapped_junc_df = mapped_junc_df[['read_name', 'aln', 'n_junctions', 'chrom', 'chain', 'donor', 'annot_donor', 'acceptor', 'annot_acceptor', 'chimeric']].sort_values(by=['read_name','aln']).reset_index(drop=True)        
gr = mapped_junc_df.groupby(['read_name','aln']).apply(lambda x: x.chimeric.any()).reset_index(name='chim_read')
del gr['aln']
mapped_junc_df = pd.merge(mapped_junc_df, gr, on='read_name').reset_index(drop=True)
mapped_junc_df.to_csv('%s/mapped_junc_df_segemehl.csv' % path_to_file, sep = '\t')
shell_call('gzip -f %s/mapped_junc_df_segemehl.csv' % path_to_file)

x = mapped_junc_df.groupby(['n_junctions','chim_read']).apply(lambda x: x.read_name.nunique()).reset_index(name='counts')
y = pd.pivot_table(x, index=['n_junctions'], columns=['chim_read'],values=['counts'], fill_value=0, aggfunc=sum, margins=True)
html_file = 'segemehl_pivot_table.html'
init_file(html_file, folder = path_to_file)
writeln_to_file(y.to_html(), html_file, folder = path_to_file)

junc_of_interest = mapped_junc_df.query('n_junctions >= 2 & chim_read == True').groupby(['read_name','aln'])
junc_csv_name = 'junc_of_interest.csv'
with open('%s/%s' % (path_to_file, junc_csv_name), 'w') as junc_csv:
    for name, group in junc_of_interest:
        print name
        print group
        junc_csv.write('\t'.join([name[0], str(name[1])]) + '\n')

                                                      