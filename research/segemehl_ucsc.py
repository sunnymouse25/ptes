# Takes json outputs from segemehl_encode.py
# Makes .BED file for UCSC genome browser

import json, os

import yaml
from interval import interval
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ptes.lib.general import init_file, writeln_to_file, shell_call
from ptes.ucsc.ucsc import list_to_dict, get_track_list
from ptes.ptes import order_interval_list

# Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="Folder with input files")
parser.add_argument("-g","--genome", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',    
                    help="Absolute path to genome file")                    
parser.add_argument("-o","--output", type=str,
                    default='bed_track_list.bed',
                    help="Output BED name, will be in ./bed/")                                       
args = parser.parse_args()

### Functions

def split_by_chimeric(lst):
    '''
    Takes list of intervals
    Cuts it by chimeric junctions (previous start > current start)
    Returns parts of the list
    '''
    strt = lst[0][0].inf
    for i, value in enumerate(lst):
        current_start = value[0].inf
        if i > 0 and current_start <= strt:   # after chimeric junction            
            read_list1 = lst[:i]
            read_list2 = lst[i:] 
            return [read_list1] + [x for x in split_by_chimeric(read_list2)]
        else:
            strt = current_start
    return [lst]
    

def get_subseq(genome, strand, chrom, s1,e1):  #cut from genome by 1-based coordinates
    if strand == '+':
        seq = genome[chrom].seq[s1-1:e1-1]
    elif strand == '-':
        seq = genome[chrom].seq[s1-1:e1-1].reverse_complement()
    return seq   
    
def splice_letters(genome, strand, chr, i1, i2):   
    '''
    STAR prints first base of donor's intron (i1) and last base of acceptor's intron (i2)
    This function returns splice site letters by their coordinates
    '''
    if strand == '+':
        donor_ss = get_subseq(genome, '+', chr, i1, (i1+2))
        acceptor_ss = get_subseq(genome, '+', chr, (i2-1), i2+1)
    elif strand == '-':
        donor_ss = get_subseq(genome, '-', chr, (i1-1), (i1+1))
        acceptor_ss = get_subseq(genome, '-', chr, i2, (i2+2))
    else: 
        return "Unknown strand"
        
    return str(donor_ss), str(acceptor_ss)    

### Main  

# Input files 
# source: https://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-from-json

read_intervals_name = args.input.rstrip('/') + '/' + 'intervals_df_segemehl.json'
print 'Reading intervals file...'
with open(read_intervals_name, 'r') as read_intervals_file:
    read_intervals = yaml.safe_load(read_intervals_file)
print 'done'
    
read_infos_name = args.input.rstrip('/') + '/' + 'infos_df_segemehl.json'  
print 'Reading infos file...'    
with open(read_infos_name, 'r') as read_infos_file:
    read_infos = yaml.safe_load(read_infos_file)
print 'done'
    
path_to_file = os.path.dirname(os.path.realpath(read_intervals_name))
if path_to_file == '':
    path_to_file = '.'

junctions_name = args.input.rstrip('/') + '/' + 'junc_of_interest.csv'
junc_of_interest = []
with open(junctions_name, 'r') as junctions_file:
    for line in junctions_file:
        junc_of_interest.append(line.strip().split('\t'))

print 'Reading genome file...'        
genome_file = args.genome
genome = SeqIO.index(genome_file, "fasta")
print 'done'
        
folder_name = '%s/bed/' % path_to_file
cmd1 = 'if [ ! -d %s ]; then mkdir %s; fi' % (folder_name, folder_name)    
shell_call(cmd1)
bed_name = args.output
coord_name = bed_name + '.coords.csv'
init_file(bed_name, folder=folder_name)   # one BED file for all tracks
init_file(coord_name, folder=folder_name)   # read_name - window in genome browser
writeln_to_file('browser full knownGene ensGene cons100way wgEncodeRegMarkH3k27ac', bed_name, folder=folder_name)
writeln_to_file('browser dense refSeqComposite pubs snp150Common wgEncodeRegDnaseClustered wgEncodeRegTfbsClusteredV3', bed_name, folder=folder_name)
writeln_to_file('browser pack gtexGene', bed_name, folder=folder_name)

for line_list in junc_of_interest:   
    key = line_list[0]  # key is mapped read_name, xi - number of current read alignment
    xi = line_list[1]
    n_junctions = line_list[2]
    annot_donors = line_list[3]
    annot_acceptors = line_list[4]    
    chrom = read_infos[key][xi][0]           
    chain = read_infos[key][xi][1]        
    tuples = read_intervals[key][xi]   # list of intervals of current read alignment    
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
    if sum(map(lambda x: x == 'GT/AG', junction_letters)) < 1:
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

        