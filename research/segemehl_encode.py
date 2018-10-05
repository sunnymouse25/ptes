# Takes segemehl.sam.nohead as input, makes junction table and HTML table for junctions counts

# Imports
import subprocess, sys, os, pickle, json, yaml, datetime
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np
from interval import interval
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ptes.lib.general import init_file, writeln_to_file, shell_call
from ptes.constants import PTES_logger
from ptes.ptes import get_read_interval, one_interval, get_interval_length, get_subseq
from ptes.ucsc.ucsc import order_interval_list, list_to_dict, get_track_list


### Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="Segemehl output file; SAM file without header")                      
parser.add_argument("-o","--output", type=str,
                    help="Output folder for results")  
parser.add_argument("-g","--genome", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',    
                    help="Absolute path to genome file")                    
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
        
def split_by_chimeric(lst):
    '''
    Takes list of intervals
    Cuts it by chimeric junctions (previous start > current start)
    Returns parts of the list
    '''
    prev_end = lst[0][0].sup
    for i, value in enumerate(lst):
        current_start = value[0].inf
        if i > 0 and current_start <= prev_end:   # after chimeric junction            
            read_list1 = lst[:i]
            read_list2 = lst[i:] 
            return [read_list1] + [x for x in split_by_chimeric(read_list2)]
        else:
            prev_end = value[0].sup
    return [lst]

def print_time():
    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M")

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
print_time()
PTES_logger.info('Reading SAM input...')
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


print_time()
PTES_logger.info('Reading SAM input... done')
    
# Exons GTF to junctions dict

PTES_logger.info('Reading GTF...')
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
print_time()
PTES_logger.info('Reading GTF... done')    
    
# Mapped junctions table 
PTES_logger.info('Creating junctions table...')           
intervals_list = []   # to remember read intervals
junc_list = []   # from mapped read intervals to list of junctions

for key in read_intervals:   # key is read_name    
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
print_time()
PTES_logger.info('Creating junctions table... done')           

x = mapped_junc_df.groupby(['n_junctions','chim_read']).apply(lambda x: x.read_name.nunique()).reset_index(name='counts')
y = pd.pivot_table(x, index=['n_junctions'], columns=['chim_read'],values=['counts'], fill_value=0, aggfunc=sum, margins=True)
html_file = 'segemehl_pivot_table.html'
init_file(html_file, folder = path_to_file)
writeln_to_file(y.to_html(), html_file, folder = path_to_file)

junc_of_interest = mapped_junc_df.query('n_junctions >= 2 & chim_read == True').sort_values(by=['annot_donor','annot_acceptor'], ascending=False).reset_index(drop=True).groupby(['read_name','aln'])
junc_csv_name = 'junc_of_interest.csv'

PTES_logger.info('Reading genome file...')        
genome_file = args.genome
genome = SeqIO.index(genome_file, "fasta")
print_time()
PTES_logger.info('Reading genome file... done')        

PTES_logger.info('Creating BED file...')
folder_name = '%s/bed/' % path_to_file
cmd1 = 'if [ ! -d %s ]; then mkdir %s; fi' % (folder_name, folder_name)    
shell_call(cmd1)
bed_name = 'bed_track_list.bed'
coord_name = bed_name + '.coords.csv'
init_file(bed_name, folder=folder_name)   # one BED file for all tracks
init_file(coord_name, folder=folder_name)   # read_name - window in genome browser
writeln_to_file('browser full knownGene ensGene cons100way wgEncodeRegMarkH3k27ac', bed_name, folder=folder_name)
writeln_to_file('browser dense refSeqComposite pubs snp150Common wgEncodeRegDnaseClustered wgEncodeRegTfbsClusteredV3', bed_name, folder=folder_name)
writeln_to_file('browser pack gtexGene', bed_name, folder=folder_name)

with open('%s/%s' % (path_to_file, junc_csv_name), 'w') as junc_csv:
    for name, group in junc_of_interest:
        key = name[0]  # key is mapped read_name, xi - number of current read alignment
        xi = name[1]
        n_junctions = group.n_junctions.iloc[0]
        annot_donors = sum(list(group.annot_donor))
        annot_acceptors = sum(list(group.annot_acceptor))
        infos = read_infos[key][str(xi)]
        chroms = set(infos[::2])
        chains = set(infos[1::2])
        chrom = chroms.pop()           
        chain = chains.pop()        
        tuples = read_intervals[key][str(xi)]   # list of intervals of current read alignment    
        tuples = sorted(tuples, key=lambda x:x[0])   # sort by xq             
        values =  [i[1] for i in tuples]   # get rid of xq
        values =  [interval[x[0][0], x[0][1]] for x in values]   # get rid of xq
        values = order_interval_list(values)   # ascending order is essential for BED lines
        junction_letters = []    
        for i in range(len(values)):
            if i > 0:
                if chain == '+':                        
                    i1 = int(values[i-1][0].sup) + 1
                    i2 = int(values[i][0].inf) - 1
                    donor_ss = get_subseq(genome, '+', chrom, i1, (i1+2))
                    acceptor_ss = get_subseq(genome, '+', chrom, (i2-1), i2+1)
                elif chain == '-':                        
                    i1 = int(values[i-1][0].inf) - 1
                    i2 = int(values[i][0].sup) + 1
                    donor_ss = get_subseq(genome, '-', chrom, (i1-1), (i1+1))
                    acceptor_ss = get_subseq(genome, '-', chrom, i2, (i2+2))
                junction = '%s/%s' % tuple(map(str, [donor_ss, acceptor_ss]))
                junction_letters.append(junction)
        track_lists = []
        windows_min = []
        windows_max = []
        if sum(map(lambda x: x == 'GT/AG', junction_letters)) < 2:
            continue        
        for i, part in enumerate(split_by_chimeric(values)):                    
            read_dict = list_to_dict(part)              
            track_list = get_track_list(chrom, chain, read_dict, name='part_%i' % (i+1))
            track_lists.append(track_list)           
            windows_min.append(int(track_list[1]))   # track_list[1] is chromStart, track_list[2] is chromEnd
            windows_max.append(int(track_list[2]))
        window = (chrom, 
                  min(windows_min)-200, 
                  max(windows_max)+200) 
        if max(windows_max) - min(windows_min) >= 1000000:   # too heavy for genome browser 
            continue
        read_name = key.replace('/','_').replace(':', '_')  
        junction_letters_str = '; '.join(junction_letters)
        description = 'letters %s; annot donors %s/%s, acceptors %s/%s' % (junction_letters_str,
                                                                            annot_donors,
                                                                            n_junctions,
                                                                            annot_acceptors,
                                                                            n_junctions)
        writeln_to_file('browser position %s:%i-%i' % window, bed_name, folder=folder_name)
        writeln_to_file(key + '\t%s:%i-%i\t' % window + description, coord_name, folder=folder_name)
        track_desc = 'track name="%s" \
        description="segemehl output; %s" \
        visibility=2 \
        itemRgb="On"' % (key, 
                        description)
        writeln_to_file(track_desc, bed_name, folder=folder_name)
        for track_list in track_lists:
            writeln_to_file('\t'.join(track_list), bed_name, folder=folder_name)
        junc_csv.write('\t'.join(map(str,[name[0], 
                                    name[1], 
                                    group.n_junctions.iloc[0],
                                    chain,
                                    sum(list(group.annot_donor)),
                                    sum(list(group.annot_acceptor)),
                                    junction_letters_str])) + '\n')   

print_time()
PTES_logger.info('Creating BED file... done')
                                                      