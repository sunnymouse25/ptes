from collections import defaultdict, OrderedDict
from interval import interval
import numpy as np


def order_cigar(row):
    """
    Returns pairs attribute - number in the order as in CIGAR: (['S1', 'M1'], [44, 32])
    """
    attr = []
    nums = []
    num = ''
    index_dict = defaultdict(lambda: 1)
    for letter in row:
        if letter in '-0123456789':
            num += letter
        if letter in 'MDINSpH':
            attr.append(letter + str(index_dict[letter]))
            index_dict[letter] += 1
            nums.append(int(num))
            num = ''
    return attr, nums


def get_read_interval(cigar, leftpos, output='dict'):
    features, lengths = order_cigar(cigar)
    if isinstance(leftpos, str):
        leftpos = int(leftpos)
    cigar_dict = OrderedDict(zip(features, lengths))  # for match_interval
    read_interval = interval()
    read_list = []
    i = 0  # counting features
    indels = []
    i_to_del = False
    for feature, length in cigar_dict.items():  # OrderedDict is important here
        if i == 0 and (feature == 'S1' or feature == 'H1'):  # soft-clip in the beginning is not included in leftpos
            k = interval[leftpos - length, leftpos - 1]
        elif 'I' in feature:  # does not consume reference
            indels.append(feature)
            i_to_del = True
            continue
        else:
            k = interval[leftpos, leftpos + length - 1]
            leftpos = leftpos + length
        read_interval = read_interval | k
        read_list.append(k)
        i += 1

    if i_to_del:  # indels do not have genomic intervals
        for x in indels:  # we need to delete them to achieve the same length of features and intervals
            del cigar_dict[x]
    read_dict = OrderedDict(zip(cigar_dict.keys(),
                                read_list))  # {'S1': interval([586.0, 657.0]), 'M1': interval([658.0, 686.0])}
    if output == 'interval':
        return read_interval
    elif output == 'dict':
        return read_dict


def read_sam_attrs(sam_name, tag_list):
    """
    Reads SAM file (no header), looks for tags in extra fields
    Example tag_list = ["NH", "HI", "XO","XS","XT"]
    Returns list ready for DataFrame:
    df_sam = pd.DataFrame(csv_list)
    """
    csv_list = []
    with open(sam_name, 'r') as sam_file:
        for line in sam_file:
            row = line.strip().split('\t')
            flag = row[1]
            if flag == 4:
                continue  # unmapped
            sam_attrs = {'read_name': row[0],
                         'flag': flag,
                         'chrom': row[2],
                         'leftpos': row[3],
                         'cigar': row[5]}
            tags = dict.fromkeys(tag_list, np.nan)
            for elm in row[11:]:
                for tag in tag_list:
                    if tag in elm:
                        tags[tag] = elm
            sam_attrs.update(tags)
            csv_list.append(sam_attrs)
    return csv_list


def one_interval(I):
    """
    Returns one large interval from set of intervals
    with extreme points of I
    """
    return interval[I[0].inf, I[-1].sup]


def get_interval_length(i):  # i is interval object of any length
    interval_length = 0
    for x in i.components:
        if len(x) > 0: interval_length += (x[0].sup - x[0].inf + 1)
    return int(interval_length)


def len_to_coord(len_list):  # from lengths of exons to the list of intervals
    intervals = []
    coord = 0  # 0-based position
    for length in len_list:
        k = interval[coord, coord + length - 1]
        coord = coord + length
        intervals.append(k)
    return intervals


def actual_exon_numbers(exons):
    """
    Takes a string from transcript name
    Returns list of actual exon numbers in transcript
    """
    exon_list = list()
    if len(exons) < 4 or len(exons) > 12:
        return [0, 0, 0, 0]
    if len(exons) == 4:
        exon_list = list(exons)
    if len(exons) == 5:
        exon_list = map(str, [7, 9, 8, 10])
    if len(exons) == 6:
        exon_list = map(str, [8, 10, 9, 11])
    if len(exons) == 7:
        exon_list = map(str, [9, 11, 10, 12])
    if len(exons) == 8:
        exon_list = [exons[0:2], exons[2:4], exons[4:6], exons[6:8]]
    if len(exons) == 9:
        exon_list = map(str, [97, 99, 98, 100])
    if len(exons) == 10:
        exon_list = map(str, [98, 100, 99, 101])
    if len(exons) == 11:
        exon_list = map(str, [99, 101, 100, 102])
    if len(exons) == 12:
        exon_list = [exons[0:3], exons[3:6], exons[6:9], exons[9:12]]
    return exon_list


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
    