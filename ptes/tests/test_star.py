

import os
import shutil
import unittest
import json
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

from ptes.star import star_SE_chimeric
from ptes.constants import TEST_DIR, PTES_logger
from ptes import ptes


PACKAGE_DIR = 'ptes'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestStar(unittest.TestCase):
    gtf_exons_name = os.path.join(INPUT_DIR,'hg19_exons_prot_coding.gtf.shuf')
    gtf_donors, gtf_acceptors = ptes.annot_junctions(gtf_exons_name=gtf_exons_name)

    chim_list = []
    with open(os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'), 'r') as chim_file:
        for line in chim_file:
            res_dict = ptes.star_line_dict(line=line)
            chim_list.append(res_dict)

    def test_parse_sam_row(self):
        """
        5 random lines from ENCFF636QII/mate1_Aligned.out.sam
        :return: dicts with attributes
        """
        with open(os.path.join(INPUT_DIR,'Aligned.out.sam.shuf')) as sam_file:
            for line in sam_file:
                row = line.strip().split('\t')
                if len(row) > 1:
                    res_dict = ptes.parse_sam_row(row=row)

    def test_sam_input(self, dump=False):
        res_dict = star_SE_chimeric.sam_input(sam_name=os.path.join(INPUT_DIR, 'Aligned.out.sam.shuf'),
                                              chim_name=os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'))
        if dump:
            with open(os.path.join(OUTPUT_DIR, 'sam_dict.json'), 'w') as sam_json_file:
                json.dump(res_dict, sam_json_file, indent=2)
        else:
            with open(os.path.join(INPUT_DIR, 'sam_dict.json'), 'r') as sam_json_file:
                res_dict_exp = json.load(sam_json_file)
            self.assertEqual(res_dict, res_dict_exp)

    def test_chim_input(self, dump=False):
        junc_dict = defaultdict(list)
        sam_dict = star_SE_chimeric.sam_input(sam_name=os.path.join(INPUT_DIR, 'Aligned.out.sam.shuf'),
                                              chim_name=os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'),
                                              )

        chim_res_list = star_SE_chimeric.chim_input(chim_name=os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'),
                                                    gtf_donors=self.gtf_donors,
                                                    gtf_acceptors=self.gtf_acceptors,
                                                    sam_dict=sam_dict,
                                                    junc_dict=junc_dict)
        if dump:
            with open(os.path.join(OUTPUT_DIR, 'chim_dict.json'), 'w') as chim_json_file:
                json.dump(chim_res_list, chim_json_file, indent=2)
        else:
            with open(os.path.join(INPUT_DIR, 'chim_dict.json'), 'r') as chim_json_file:
                res_list_exp = json.load(chim_json_file)
            self.assertEqual(chim_res_list, res_list_exp)



    def test_reads_to_junctions(self, dump=False):
        res_dict = star_SE_chimeric.sam_input(sam_name=os.path.join(INPUT_DIR, 'Aligned.out.sam.shuf'),
                                              chim_name=os.path.join(INPUT_DIR,'Chimeric.out.junction.shuf'),
                                              )

        reads_list = star_SE_chimeric.chim_input(chim_name=os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'),
                                                 gtf_donors=self.gtf_donors,
                                                 gtf_acceptors=self.gtf_acceptors,
                                                 sam_dict=res_dict)


        reads_df = pd.DataFrame(reads_list)
        reads_df['id'] = 'tag'
        junc_df = star_SE_chimeric.reads_to_junctions(reads_df)
        for index, row in junc_df.head(n=3).iterrows():
            print(list(row.index))
            print(row.values)
        if dump:
            junc_df.to_csv(os.path.join(OUTPUT_DIR, 'junctions.csv'), sep='\t')
        else:
            exp_junc_df = pd.read_csv(os.path.join(INPUT_DIR, 'junctions.csv'), sep='\t')
            self.assertEqual(exp_junc_df.to_dict('records'),
                             junc_df.reset_index(drop=False).to_dict('records'))

    def test_junc_dict(self):
        junc_dict_test = defaultdict(list)

        res_dict = star_SE_chimeric.sam_input(sam_name=os.path.join(INPUT_DIR, 'Aligned.out.sam.shuf'),
                                              chim_name=os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'),
                                              )
        reads_list = star_SE_chimeric.chim_input(chim_name=os.path.join(INPUT_DIR, 'Chimeric.out.junction.shuf'),
                                                 gtf_donors=self.gtf_donors,
                                                 gtf_acceptors=self.gtf_acceptors,
                                                 sam_dict=res_dict,
                                                 junc_dict=junc_dict_test)

        with open(os.path.join(OUTPUT_DIR, 'junc_dict.json'), 'w') as junc_json:
            json.dump({str(k): v for k, v in junc_dict_test.items()}, junc_json, indent=2)

        data = json.load(open(os.path.join(OUTPUT_DIR, 'junc_dict.json')), object_pairs_hook=OrderedDict)

        with open(os.path.join(OUTPUT_DIR, 'junc_dict_loaded.json'), 'w') as junc_json:
            json.dump(data, junc_json, indent=2)
        





