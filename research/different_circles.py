# Run at /home/sunnymouse/projects/PTES/panhandles

import pandas as pd
import numpy as np


all_circles = pd.read_csv('all_circles/bed/all_Chimeric.single.bed.sorted', sep='\t',header=None)
cols = ['chrom','start','end','name','score','chain','th_start','th_end','color','blocks','lens','starts']
all_circles.columns = cols

junc = pd.read_csv('../ENCODE/chim_types.csv', sep='\t', index_col = 0)

inside_list = []
inside_list_unique = []
outside_list = []
outside_list_unique = []
unknown_list = []
unknown_list_unique = []

inside_unique_set = set()
outside_unique_set = set()
unknown_unique_set = set()

for i, row in all_circles.iterrows():
    chrom = row['chrom']
    name = row['name']
    name_list = name.split('_')
    donor_ss = name_list[1]
    acceptor_ss = name_list[2]
    q = junc.query('chrom == "%s" & donor == %s & acceptor == %s' % (chrom, donor_ss, acceptor_ss))
    if q.shape[0] > 0:
        if q['outside'].any() > 0:
            outside_list.append(row)
            if (chrom, donor_ss, acceptor_ss) not in outside_unique_set:
                outside_unique_set.add((chrom, donor_ss, acceptor_ss))
                outside_list_unique.append(row)
        if q['inside'].any() > 0:
            inside_list.append(row)
            if (chrom, donor_ss, acceptor_ss) not in inside_unique_set:
                inside_unique_set.add((chrom, donor_ss, acceptor_ss))
                inside_list_unique.append(row)

    else:
        unknown_list.append(row)
        if (chrom, donor_ss, acceptor_ss) not in unknown_unique_set:
            unknown_unique_set.add((chrom, donor_ss, acceptor_ss))
            unknown_list_unique.append(row)


inside_df = pd.DataFrame(inside_list)
inside_df.to_csv('all_Chimeric.inside.single.bed', sep='\t', index=False, header=False)

inside_df_unique = pd.DataFrame(inside_list_unique)
inside_df_unique.to_csv('all_Chimeric.inside.unique.single.bed', sep='\t', index=False, header=False)

outside_df = pd.DataFrame(outside_list)
outside_df.to_csv('all_Chimeric.outside.single.bed', sep='\t', index=False, header=False)

outside_df_unique = pd.DataFrame(outside_list_unique)
outside_df_unique.to_csv('all_Chimeric.outside.unique.single.bed', sep='\t', index=False, header=False)

unknown_df = pd.DataFrame(unknown_list)
unknown_df.to_csv('all_Chimeric.unknown.single.bed', sep='\t', index=False, header=False)

unknown_df_unique = pd.DataFrame(unknown_list_unique)
unknown_df_unique.to_csv('all_Chimeric.unknown.unique.single.bed', sep='\t', index=False, header=False)

