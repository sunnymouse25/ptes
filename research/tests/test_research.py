import os
import unittest

import ptes.ptes
from research.constants import TEST_DIR

PACKAGE_DIR = 'research'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestStarEncode(unittest.TestCase):
    def test_annot_junctions(self):
        filename = './test_data/hg19_exons_prot_coding.gtf.shuf'
        donors, acceptors = ptes.ptes.annot_junctions(gtf_exons_name=filename)
        chrom = 'chr16'
        real_donors = [11870307, 12093216]
        real_acceptors = [11870215, 12093153]
        for real_donor in real_donors:
            self.assertIn(real_donor, donors[chrom])
        for real_acceptor in real_acceptors:
            self.assertIn(real_acceptor, acceptors[chrom])

if __name__ == "__main__":
    unittest.main()