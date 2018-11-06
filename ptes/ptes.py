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
    if not isinstance(leftpos, int):
        leftpos = int(leftpos)
    cigar_list = zip(features, lengths)  # for match_interval
    read_interval = interval()
    read_list = []    
    features_list = []
    for feature, length in cigar_list:  
        if 'M' in feature or 'D' in feature or 'N' in feature or 'p' in feature:
            k = interval[leftpos, leftpos + length - 1]
            leftpos = leftpos + length
            read_interval = read_interval | k
            read_list.append(k) 
            features_list.append(feature)

    read_dict = OrderedDict(zip(features_list,
                                read_list))  # {'S1': interval([586.0, 657.0]), 'M1': interval([658.0, 686.0])}
    if output == 'interval':
        return read_interval
    elif output == 'dict':
        return read_dict
    elif output == 'list':
        return read_list


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

def split_by_p(read_dict):
    '''
    Takes read_dict (OrderedDict),
    returns list of dicts: before p and after p
    '''
    items = read_dict.items()
    for i, item in enumerate(items):
        if 'p' in item[0]:
            return [OrderedDict(items[:i]), OrderedDict(items[(i+1):])]
    return [read_dict]


def annot_junctions(gtf_exons_name):
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
                gtf_donors[chrom].add(end + 1)
                gtf_acceptors[chrom].add(strt - 1)
            elif chain == '-':
                gtf_donors[chrom].add(strt - 1)
                gtf_acceptors[chrom].add(end + 1)
    return gtf_donors, gtf_acceptors


def dict_to_interval(read_dict, put_n = True, output='interval'):
    '''
    Takes read_dict (OrderedDict),
    returns interval of intervals for all features that consume reference
    '''
    output_interval = interval()
    for item in read_dict.items():
        feature = item[0]
        if put_n:
            if 'M' in feature or 'D' in feature or 'N' in feature:
                output_interval = output_interval | item[1]
        else:
            if 'M' in feature or 'D' in feature:
                output_interval = output_interval | item[1]
    if output == 'interval':
        return output_interval
    if output == 'list':
        return [x for x in output_interval.components]


def mate_intersection(interval1, interval2):
    intersection = one_interval(interval1) & one_interval(interval2)
    if intersection == interval():  # zero intersection
        return 'outside'
    else:
        return 'inside'


def return_mates(cigar1, coord1, cigar2, coord2, chain):
    mate1 = None
    mate2 = None
    chim_part1 = get_read_interval(cigar1, coord1)  # not mates, chimeric parts!
    chim_part2 = get_read_interval(cigar2, coord2)
    if 'p' in cigar1:
        splits = split_by_p(chim_part1)
        if chain == '+':
            mate2 = one_interval(dict_to_interval(splits[0]))
            mate_intervals = dict_to_interval(splits[1]) | dict_to_interval(chim_part2)
            mate1 = one_interval(mate_intervals)
        if chain == '-':
            mate_intervals = dict_to_interval(splits[0]) | dict_to_interval(chim_part2)
            mate1 = one_interval(mate_intervals)
            mate2 = one_interval(dict_to_interval(splits[1]))
    elif 'p' in cigar2:
        splits = split_by_p(chim_part2)
        if chain == '+':
            mate_intervals = dict_to_interval(chim_part1) | dict_to_interval(splits[0])
            mate1 = one_interval(mate_intervals)
            mate2 = one_interval(dict_to_interval(splits[1]))
        if chain == '-':
            mate_intervals = dict_to_interval(chim_part1) | dict_to_interval(splits[1])
            mate1 = one_interval(mate_intervals)
            mate2 = one_interval(dict_to_interval(splits[0]))
    return mate1, mate2


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


def interval_to_string(i):
    """
    Useful for writing intervals to file
    :param i: i is interval object with length 1
    :return: string interpretation of given interval
    """
    return str([str(int(i[0].inf)), str(int(i[0].sup))])


def sort_by_xq(tuples, chain):
    """
    Takes value from read_intervals dictionary
    Example for value: [(0,A), (1,B), (2,C)],
    where A,B,C are M-intervals of non-chim read
    Chain + 5'-3' A-B-C
    Chain - 5'-3' C-B-A
    Makes list of M-intervals in the order as they are in the read,
    but from 5' to 3', as in STAR Chimeric.out.junction
    :param tuples: example for chain -
    [(0, interval([74603869.0, 74603880.0])),
     (1, interval([74600055.0, 74600075.0])),]
    :param chain: if chain == '-' than their order is from 3' to 5'
    :return: list of M-intervals in the ascending order 5'-3'
    """
    tuples = sorted(tuples, key=lambda x: x[0])  # sort by xq
    values = [i[1] for i in tuples]  # get rid of xq
    if chain == '-':
        values = values[::-1]
    return values


def get_junctions(chrom, chain, values, gtf_donors, gtf_acceptors):
    """
    From list of segemehl M-intervals to junctions
    :param chrom: chromosome
    :param chain: chain relative to genome
    :param values: sorted by xq list of M-intervals ordered from 5' to 3'
    :param gtf_donors: dictionary with chromosomes as keys, sets of ints as values
    :param gtf_acceptors: same as gtf_donors, output of annot_junction function
    :return: list of dictionaries with coordinates of junctions, ready to make pandas dataframe
    """
    n_j = len(values) - 1
    junc_list = []  # from mapped read intervals to the list of junctions
    for i in range(1, len(values)):
        donor_ss = np.nan
        acceptor_ss = np.nan
        chimeric = False
        if chain == '+':
            donor_ss = int(values[i-1][0].sup) + 1
            acceptor_ss = int(values[i][0].inf) - 1
            if donor_ss > acceptor_ss:
                chimeric = True
        elif chain == '-':
            donor_ss = int(values[i][0].inf) - 1
            acceptor_ss = int(values[i-1][0].sup) + 1
            if donor_ss < acceptor_ss:
                chimeric = True
        donors = gtf_donors[chrom]
        acceptors = gtf_acceptors[chrom]
        annot_donor = 1 if donor_ss in donors else 0
        annot_acceptor = 1 if acceptor_ss in acceptors else 0
        junc_list.append({'n_junctions': n_j,
                'chrom': chrom,
                'chain': chain,
                'donor': str(donor_ss),
                'annot_donor': annot_donor,
                'acceptor': str(acceptor_ss),
                'annot_acceptor': annot_acceptor,
                'chimeric': chimeric})
    return junc_list