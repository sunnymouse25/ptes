from collections import OrderedDict
import os

from interval import interval

from ptes.ptes import get_interval_length, dict_to_interval
from ptes.lib.general import init_file, writeln_to_file, shell_call, make_dir


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
    Returns .bed line with all tracks in the list, from 1-based to 0-based
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


def make_bed_folder(prefix, path_to_file):
    """
    Initiates 3 files essential for Genome Browser:
    :param path_to_file: folder where ./bed subfolder will be
    :param prefix: prefix for names of files, i.e. sample tag
    :return: ./bed subfolder,
    BED file for track lines (.bed),
    table with windows to copy-paste (.coords.csv),
    track file for GB (.track)
    """
    bed_name = '%s.bed' % prefix  # only track lines
    coord_name = '%s.coords.csv' % prefix  # table with windows to paste into GB and with descriptions
    info_name = '%s.track' % prefix  # file to submit to GB
    folder_name = '%s/bed/' % path_to_file
    make_dir(folder_name)

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
    prefix, bed_name.replace('.bed', '.bb')), info_name, folder=folder_name)
    return folder_name, bed_name, coord_name


def to_bigbed(bed_name, folder_name):
    """
    Runs bedToBigBed script to convert bed to bigBed,
    bedToBigBed must be in $PATH
    :param bed_name: Name of BED file to be converted
    :param folder_name: Folder of BED file
    :return: sorted bed and bigBed files in the same folder
    """
    if folder_name[-1] != '/':
        folder_name = folder_name + '/'
    real_path = os.path.dirname(os.path.realpath(folder_name+bed_name))
    cmd1 = 'sort -k1,1 -k2,2n %s/%s > %s/%s.sorted' % (real_path, bed_name, real_path, bed_name)
    cmd2 = 'bedToBigBed \
            %s/%s.sorted \
            http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes \
            %s/%s' % (real_path, bed_name, real_path, bed_name.rpartition('.')[0]+'.bb')
    cmds = [cmd1,cmd2]
    for cmd in cmds:
        shell_call(cmd)


def single_track(read_dict_list, kwargs):
    """
    Makes single BED line for the list of read_dicts
    :param read_dict_list: list of read_dicts (get_read_interval outputs)
    :param kwargs: kwargs for get_track_list function
    :return: track list ready for .bed output
    """
    interval_list = map(lambda x: dict_to_interval(x, put_n=False), read_dict_list)
    single_interval = interval()
    for part in interval_list:
        single_interval = single_interval | part
    single_interval_list = [y for y in single_interval.components]
    single_track = get_track_list(
                read_dict=list_to_dict(single_interval_list), **kwargs)  # for checking in GB that intervals are same
    return single_track