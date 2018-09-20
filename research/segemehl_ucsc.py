# Takes json outputs from segemehl_encode.py
# Makes .BED file for UCSC genome browser

import json, os

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

### Main  

# Input files 
read_intervals_name = args.input.rstrip('/') + '/' + 'intervals_df_segemehl.json'
with open(read_intervals_name, 'r') as read_intervals_file:
    read_intervals = json.load(read_intervals_file)
    
read_infos_name = args.input.rstrip('/') + '/' + 'infos_df_segemehl.json'  
with open(read_infos_name, 'r') as read_infos_file:
    read_infos = json.load(read_infos_file)
    
path_to_file = os.path.dirname(os.path.realpath(read_intervals_name))
if path_to_file == '':
    path_to_file = '.'

junctions_name = args.input.rstrip('/') + '/' + 'junc_of_interest.csv'
with open(junctions_name, 'r') as junctions_file:
    junc_of_interest = set(junctions_file.read().split('\n'))
    
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


for key in read_intervals.keys():   # key is mapped read_name    
    if key not in junc_of_interest:
        continue            
    xi = 0
    while read_intervals[key].get(str(xi), None):   # xi - number of current read alignment                
        chrom = read_infos[key][str(xi)][0]           
        chain = read_infos[key][str(xi)][1]        
        tuples = read_intervals[key][str(xi)]   # list of intervals of current read alignment    
        tuples = sorted(tuples, key=lambda x:x[0])   # sort by xq        
        values =  [i[1] for i in tuples]   # get rid of xq
        values =  [interval[x[0][0], x[0][1]] for x in values]   # get rid of xq
        values = order_interval_list(values)   # ascending order is essential for BED lines
        chimeric = False    # by default no chimeric junctions assumed
        start = values[0][0].inf
        for i, value in enumerate(values):
            current_start = value[0].inf
            if current_start < start:   # after chimeric junction
                read_list1 = values[:i]
                read_list2 = values[i:]
                chimeric = True
                break
            else:
                start = current_start
        xi += 1
        if chimeric:
            read_dict1 = list_to_dict(read_list1)
            read_dict2 = list_to_dict(read_list2)                
            if chrom.startswith('chr'):                                                 
                track_list1 = get_track_list(chrom, chain, read_dict1, name='chim_donor')
                track_list2 = get_track_list(chrom, chain, read_dict2, name='chim_acceptor')            
                junction_letters = 'NA'                        
                window = (chrom, 
                      min(int(track_list1[1]), int(track_list2[1]))-200, 
                      max(int(track_list1[2]), int(track_list2[2]))+200)    # track_list[1] is chromStart, track_list[2] is chromEnd

                read_name = key.replace('/','_').replace(':', '_')                            
                writeln_to_file('browser position %s:%i-%i' % window, bed_name, folder=folder_name)
                writeln_to_file(key + '\t%s:%i-%i' % window, coord_name, folder=folder_name)
                track_desc = 'track name="%s" description="segemehl output" visibility=2 itemRgb="On"' % key       
                writeln_to_file(track_desc, bed_name, folder=folder_name)
                writeln_to_file('\t'.join(track_list1), bed_name, folder=folder_name)
                writeln_to_file('\t'.join(track_list2), bed_name, folder=folder_name)