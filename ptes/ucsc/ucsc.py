from collections import OrderedDict, defaultdict

from ptes.ptes import get_interval_length


def list_to_dict(lst):
    """
    From list of intervals makes dict of M's
    {'M1': [interval([70799946.0, 70799962.0], 'M2': interval([70801644.0, 70801667.0])}
    """
    numbers = range(1, len(lst)+1)
    names = map(lambda x: 'M' + str(x), numbers)
    dct = OrderedDict(zip(names,lst))
    return dct


def order_interval_list(values):
    """
    For list of intervals returns list in increasing order:
    inverts for '-' chain
    """
    start = values[0][0].inf
    last_start = values[-1][0].inf  
    if start > last_start:
        values = values[::-1]
    return values


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
