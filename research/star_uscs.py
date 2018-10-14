#Input: row(s) from Chimeric.out.junction
#Output: file(s) for viewing in UCSC genome browser, one per row
#Source for BED format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#Usage: python star_ucsc.py -i list -o output_prefix -d description (PTES type)


###IMPORT
import os
from collections import defaultdict, OrderedDict

from interval import interval
import pandas as pd

from ptes.lib.general import init_file, writeln_to_file, shell_call
from ptes.ptes import get_read_interval, one_interval, get_interval_length, split_by_p
from ptes.ucsc.ucsc import order_interval_list, list_to_dict, get_track_list


### Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="Name of input file")
                    
parser.add_argument("-o","--output", type=str,
                    default = 'output.bed',
                    help="Prefix of output file, i.e. ENCODE id")
                    
parser.add_argument("-d","--description", type=str,                    
                    help="Description for tracks in BED (PTES type)")                    
args = parser.parse_args()                    


### Functions
  
def split_cigar(cigar, leftpos): 
    '''
    Returns cigar and leftpos for the 2nd mate
    '''      
    num = ''    
    new_leftpos = leftpos
    new_cigar = cigar.split('p')[1] 
    for letter in cigar:
        if letter in '-0123456789':
            num += letter
        elif letter in 'MDN':  # consume reference   
            new_leftpos +=(int(num))
            num = ''
        elif letter == 'p':
            new_leftpos +=(int(num))
            return new_cigar, new_leftpos
        elif letter in 'SI':   # first S should not be included, last S is after p
            num = '0'                       
    

    
### Main

# input files
input_name = args.input    
output_prefix = args.output
desc = args.description

dirname = os.path.dirname(os.path.realpath(input_name))
if dirname == '':
    dirname = '.'

folder_name = '%s/bed/' % dirname
cmd = 'if [ ! -d %s ]; then mkdir %s; fi' % (folder_name, folder_name)    
shell_call(cmd)
    
bed_name = '%s.bed' % output_prefix
coord_name = '%s.coords.csv' % output_prefix
info_name = '%s.track' % output_prefix
init_file(bed_name, folder = folder_name)
init_file(coord_name, folder = folder_name)
init_file(info_name, folder = folder_name)

writeln_to_file('browser full knownGene ensGene cons100way wgEncodeRegMarkH3k27ac', info_name, folder = folder_name)
writeln_to_file('browser dense refSeqComposite pubs snp150Common wgEncodeRegDnaseClustered wgEncodeRegTfbsClusteredV3', 
                info_name, 
                folder = folder_name)
writeln_to_file('browser pack gtexGene', info_name, folder = folder_name)
writeln_to_file('track type=bigBed \
                name="%s" \
                description="bigBed" \
                visibility=2 \
                itemRgb="On" \
bigDataUrl=https://github.com/sunnymouse25/ptes/blob/dev/research/bed/%s?raw=true' % (desc, bed_name.replace('.bed', '.bb')), info_name, folder = folder_name)

# start processing input file
input_df = pd.read_csv(input_name, sep='\t', header=None, dtype=str)
groups = input_df.groupby(by = [0,1,4])  #chr - junction donor coord - junction acceptor coord ; no chain
for i, group in enumerate(groups):
        print i
        print group
        line_list = group[1].iloc[0].tolist()
        num_reads = group[1].shape[0]      
        chrom = line_list[0]
        intron1 = int(line_list[1])
        chain = line_list[2]
        intron2 = int(line_list[4])
        junction_type = line_list[6]   #junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC 
        if junction_type == '1':
            junction_letters = 'GT/AG'
        elif junction_type == '2':
            junction_letters = 'CT/AC'            
        else:
            junction_letters = 'unknown'
        read_name = line_list[9]
        coord1 = int(line_list[10])
        cigar1 = line_list[11]
        coord2 = int(line_list[12])
        cigar2 = line_list[13]
        track_lists = []
        windows_min = []
        windows_max = []
        if 'p' in cigar1:
            mate1 = get_read_interval(cigar1.split('p')[0].rstrip('-0123456789'), coord1)
            mate2 = get_read_interval(*split_cigar(cigar1, coord1))   
            chim_part2 = get_read_interval(cigar2, coord2)            
            bed1 = get_track_list(chrom, chain, mate1, name='mate1', color='r')
            bed2 = get_track_list(chrom, chain, mate2, name='mate2', color='r')
            bed3 = get_track_list(chrom, chain, chim_part2, name='chim', color='r')
            track_lists = [bed1, bed2, bed3]
        elif 'p' in cigar2:
            mate1 = get_read_interval(cigar2.split('p')[0].rstrip('-0123456789'), coord2)
            mate2 = get_read_interval(*split_cigar(cigar2, coord2))   
            chim_part2 = get_read_interval(cigar1, coord1)
            bed1 = get_track_list(chrom, chain, mate1, name='mate1', color='r')
            bed2 = get_track_list(chrom, chain, mate2, name='mate2', color='r')
            bed3 = get_track_list(chrom, chain, chim_part2, name='chim', color='r')
            track_lists = [bed1, bed2, bed3]
        else:   #single-read mode
            chim_part1 = get_read_interval(cigar1, coord1)   # not mates, chimeric parts!
            chim_part2 = get_read_interval(cigar2, coord2)
            bed1 = get_track_list(chrom, chain, chim_part1, name='chim_part1', color='r')
            bed2 = get_track_list(chrom, chain, chim_part2, name='chim_part2', color='b')
            track_lists = [bed1, bed2]
            
        for track_list in track_lists:
            windows_min.append(int(track_list[1]))   # track_list[1] is chromStart, track_list[2] is chromEnd
            windows_max.append(int(track_list[2]))  
            writeln_to_file('\t'.join(track_list), bed_name, folder = folder_name)
        window = (chrom, 
                  min(windows_min)-200, 
                  max(windows_max)+200)   # bed1[2] is chromEnd        
                                      
        description = "%s; %i read(s) found; chim junction %s" % (desc, num_reads, junction_letters)
        writeln_to_file('browser position %s:%i-%i\t' % window + description, coord_name, folder = folder_name) 
        if i == 0:
            writeln_to_file('browser position %s:%i-%i\t' % window, info_name, folder = folder_name) 

        
        
        
        
        
        
        