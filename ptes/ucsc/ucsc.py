from collections import OrderedDict
import os

from ptes.ptes import get_interval_length
from ptes.lib.general import init_file, writeln_to_file, shell_call


def list_to_dict(lst):
    """
    From list of intervals makes dict of M's
    {'M1': [interval([70799946.0, 70799962.0], 'M2': interval([70801644.0, 70801667.0])}
    """
    numbers = range(1, len(lst)+1)
    names = map(lambda x: 'M' + str(x), numbers)
    dct = OrderedDict(zip(names,lst))
    return dct


def get_track_list(chrom, chain, read_dict, name = 'sim_read', color = 'r'):
    """
    Gets chrom, chain, dictionary of read intervals instead of cigar and leftpos;
    Example for read_dict: {'S1': interval([586.0, 657.0]), 'M1': interval([658.0, 686.0])} 
    Name may be CIGAR string or read_name;
    Color: red by default, may be 'b', 'g' or in RGB code like '255,0,255'
    Returns .bed line with all tracks in the list
    """
    track_list = [chrom]    
    blockSizes = []
    chromStarts = []
    i = 0
    for key, value in read_dict.items():         
        if i == 0:
            chromStart = int(value[0].inf)-1
        if key == 'M1':
            thickStart = int(value[0].inf)-1
        if 'M' in key:
            thickEnd = int(value[0].sup)
            block_len = get_interval_length(value)
            blockSizes.append(block_len)
            chromStarts.append(int(value[0].inf)-1-chromStart)            
        if 'S' in key:
            block_len = get_interval_length(value)
            blockSizes.append(block_len)
            chromStarts.append(int(value[0].inf)-1-chromStart)
        if i == len(read_dict)-1:
            chromEnd = int(value[0].sup)
        i+=1    
    track_list.append(chromStart) # chromStart, counts from 0
    track_list.append(chromEnd)   # chromEnd, not included
    track_list.append(name)      # name
    track_list.append('0')        # score
    track_list.append(chain)      # chain
    track_list.append(thickStart) # thickStart
    track_list.append(thickEnd)   # thickEnd
    if color == 'r':
        track_list.append('255,0,0')  #color
    elif color == 'g':
        track_list.append('0,255,0')  
    elif color == 'b':
        track_list.append('0,0,255') 
    else:
        track_list.append(color)        # fuchsia is '255,0,255'
    track_list.append(len(blockSizes))        # blockCount is length of blockSizes
    track_list.append(','.join(map(str,blockSizes)))         # blockSizes
    track_list.append(','.join(map(str,chromStarts)))        # chromStarts
    return map(str, track_list)


def make_bed_folder(folder_name, bed_name, coord_name, info_name, data_desc):
    """
    Initiates 3 files essential for Genome Browser
    :param folder_name: ./bed subfolder
    :param bed_name: only track lines
    :param coord_name: table with windows to paste into GB and with descriptions
    :param info_name: file to submit to GB
    :param data_desc: description for the whole dataset
    :return: ./bed subfolder, BED file for track lines, table with windows to copy-paste, track file for GB
    """

    cmd1 = 'if [ ! -d %s ]; then mkdir %s; fi' % (folder_name, folder_name)
    shell_call(cmd1)
    init_file(bed_name, folder=folder_name)
    init_file(coord_name, folder=folder_name)
    init_file(info_name, folder=folder_name)

    writeln_to_file('browser full knownGene ensGene cons100way wgEncodeRegMarkH3k27ac rmsk', info_name, folder=folder_name)
    writeln_to_file(
        'browser dense refSeqComposite pubs snp150Common wgEncodeRegDnaseClustered wgEncodeRegTfbsClusteredV3',
        info_name,
        folder=folder_name
    )
    writeln_to_file('browser pack gtexGene', info_name, folder=folder_name)
    writeln_to_file('track type=bigBed \
                    name="%s" \
                    description="bigBed" \
                    visibility=2 \
                    itemRgb="On" \
    bigDataUrl=https://github.com/sunnymouse25/ptes/blob/dev/research/bed/%s?raw=true' % (
    data_desc, bed_name.replace('.bed', '.bb')), info_name, folder=folder_name)

def to_bigbed(bed_name, folder_name):
    if folder_name[-1] != '/':
        folder_name = folder_name + '/'
    real_path = os.path.dirname(os.path.realpath(folder_name+bed_name))
    cmd1 = 'sort -k1,1 -k2,2n %s/%s > %s/input.bed' % (real_path, bed_name, real_path)
    cmd2 = 'bedToBigBed \
            %s/input.bed \
            http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes \
            %s/%s' % (real_path, real_path, bed_name.rstrip('.bed')+'.bb')
    cmds = [cmd1,cmd2]
    for cmd in cmds:
        shell_call(cmd)
