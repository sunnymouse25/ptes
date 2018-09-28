# Takes json outputs from segemehl_encode.py
# Makes .BED file for UCSC genome browser

import json, os, yaml

from interval import interval

from ptes.lib.general import init_file, writeln_to_file, shell_call
from ptes.ucsc.ucsc import order_interval_list, list_to_dict, get_track_list

# Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="Folder with input files")                    
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

### Main  

# Input files 
# source: https://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-from-json

read_intervals_name = args.input.rstrip('/') + '/' + 'intervals_df_segemehl.json'
with open(read_intervals_name, 'r') as read_intervals_file:
    read_intervals = yaml.safe_load(read_intervals_file)
    
read_infos_name = args.input.rstrip('/') + '/' + 'infos_df_segemehl.json'  
with open(read_infos_name, 'r') as read_infos_file:
    read_infos = yaml.safe_load(read_infos_file)
    
path_to_file = os.path.dirname(os.path.realpath(read_intervals_name))
if path_to_file == '':
    path_to_file = '.'

junctions_name = args.input.rstrip('/') + '/' + 'junc_of_interest.csv'
junc_of_interest = []
with open(junctions_name, 'r') as junctions_file:
    for line in junctions_file:
        junc_of_interest.append(line.strip().split('\t'))
    
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
i = 0

for line_list in junc_of_interest:   
    key = line_list[0]  # key is mapped read_name, xi - number of current read alignment
    xi = line_list[1]
    annot_donors = line_list[2]
    annot_acceptors = line_list[3]    
    chrom = read_infos[key][xi][0]           
    chain = read_infos[key][xi][1]        
    tuples = read_intervals[key][xi]   # list of intervals of current read alignment    
    tuples = sorted(tuples, key=lambda x:x[0])   # sort by xq        
    values =  [i[1] for i in tuples]   # get rid of xq
    values =  [interval[x[0][0], x[0][1]] for x in values]   # get rid of xq
    values = order_interval_list(values)   # ascending order is essential for BED lines     
    track_lists = []
    windows_min = []
    windows_max = []
    for i, part in enumerate(split_by_chimeric(values)):        
        read_dict = list_to_dict(part)              
        track_list = get_track_list(chrom, chain, read_dict, name='part_%i' % (i+1))
        track_lists.append(track_list)
        junction_letters = 'NA'   # to do: add genome slicer
        windows_min.append(int(track_list[1]))   # track_list[1] is chromStart, track_list[2] is chromEnd
        windows_max.append(int(track_list[2]))
    window = (chrom, 
              min(windows_min)-200, 
              max(windows_max)+200) 
    if max(windows_max) - min(windows_min) >= 1000000:   # too heavy for genome browser 
        continue
    read_name = key.replace('/','_').replace(':', '_')                            
    writeln_to_file('browser position %s:%i-%i' % window, bed_name, folder=folder_name)
    writeln_to_file(key + '\t%s:%i-%i' % window, coord_name, folder=folder_name)
    track_desc = 'track name="%s" description="segemehl output" visibility=2 itemRgb="On"' % key       
    writeln_to_file(track_desc, bed_name, folder=folder_name)
    for track_list in track_lists:
        writeln_to_file('\t'.join(track_list), bed_name, folder=folder_name)

        