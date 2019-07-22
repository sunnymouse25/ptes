# from ipsa output .ssj.tsv to .BED6 file
# from junctions [a,b] to intervals: makes [a-10, b+10]
# /home/sunnymouse/projects/PTES/ENCODE/bed/ipsa/all_polya_plu.bed

import os

ipsa_name = 'all_polya_plus.ssj.tsv'
#ipsa_folder = '/home/dp/ipsa/hg19/ENCODE/'
bed_name = ipsa_name.split('.')[0] + '.all_reads.bed'

with open(ipsa_name, 'r') as ipsa_file, open(bed_name, 'w') as bed_file:
    for line in ipsa_file:
        line_list = line.strip().split('\t')
        coord = line_list[0]
        coord_list = coord.split('_')
        n_reads = int(line_list[1])
        score = n_reads if n_reads < 1000 else 1000
        for i in range(n_reads):
            bed_list = list(map(str,
                                [coord_list[0],
                                 int(coord_list[1]) - 10,
                                 int(coord_list[2]) + 10,
                                 '.',
                                 score,
                                 coord_list[3]]))
            bed_file.write('\t'.join(bed_list) + '\n')
