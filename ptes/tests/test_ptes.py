import os
import shutil
import unittest
from collections import OrderedDict

from interval import interval

from ptes import ptes
from ptes.constants import TEST_DIR, PTES_logger


PACKAGE_DIR = 'ptes'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestPtes(unittest.TestCase):

    input_read_dicts = [
        OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                     ('M2', interval([152.0, 171.0])), ('p1', interval([172.0, 1171.0])),
                     ('M3', interval([1172.0, 1247.0]))]),
        OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                     ('M2', interval([152.0, 171.0])), ('p1', interval([172.0, 1171.0])),
                     ('S1', interval([1172.0, 1247.0])), ('M3', interval([1248.0, 1260.0]))]),
        OrderedDict([('S1', interval([96.0, 99.0])),
                     ('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                     ('M2', interval([152.0, 171.0])), ('p1', interval([172.0, 1171.0])),
                     ('S2', interval([1172.0, 1247.0])), ('M3', interval([1248.0, 1260.0])),
                     ('S3', interval([1261.0, 1265.0]))]),
    ]

    def test_get_read_interval(self):
        cigars = [
            '20M30S',
            '30S20M',
            '30S20M30S',
            '20M2I20M',
            '20M2D20M',
            '20M32N20M',
            '20S20M32N20M',
            '20S20M32N20M1000p76M',
        ]
        leftpos = 100
        exp_read_dicts = [
            OrderedDict([('M1', interval([100.0, 119.0])),]),
            OrderedDict([('M1', interval([100.0, 119.0])),]),
            OrderedDict([('M1', interval([100.0, 119.0])),]),
            OrderedDict([('M1', interval([100.0, 119.0])), ('M2', interval([120.0, 139.0]))]),
            OrderedDict([('M1', interval([100.0, 119.0])), ('D1', interval([120.0, 121.0])),
                         ('M2', interval([122.0, 141.0]))]),
            OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                         ('M2', interval([152.0, 171.0]))]),
            OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                         ('M2', interval([152.0, 171.0]))]),
            OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                         ('M2', interval([152.0, 171.0])), ('p1', interval([172.0, 1171.0])),
                         ('M3', interval([1172.0, 1247.0]))]),
        ]
        res_read_dicts = [ptes.get_read_interval(cigar=cigar, leftpos=leftpos, output='dict') for cigar in cigars]
        for i in range(len(exp_read_dicts)):
            self.assertEqual(exp_read_dicts[i], res_read_dicts[i])

    def test_split_by_p(self):
        read_dicts = [
            OrderedDict(
                [('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])), ('M2', interval([152.0, 171.0])),
                 ('p1', interval([172.0, 1171.0])), ('M3', interval([1172.0, 1247.0]))]),
        ]
        results = [ptes.split_by_p(read_dict=read_dict) for read_dict in self.input_read_dicts]
        exp_results = [
            [
            OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                         ('M2', interval([152.0, 171.0]))]),
            OrderedDict([('M3', interval([1172.0, 1247.0]))])
            ],
            [
                OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                             ('M2', interval([152.0, 171.0]))]),
                OrderedDict([('M3', interval([1248.0, 1260.0]))])
            ],
            [
                OrderedDict([('M1', interval([100.0, 119.0])), ('N1', interval([120.0, 151.0])),
                             ('M2', interval([152.0, 171.0]))]),
                OrderedDict([('M3', interval([1248.0, 1260.0]))])
            ],
        ]
        for i in range(len(read_dicts)):
            self.assertEqual(exp_results[i], results[i])

    def test_annot_junctions(self):
        filename = './test_data/ptes/hg19_exons_prot_coding.gtf.shuf'
        donors, acceptors = ptes.annot_junctions(gtf_exons_name=filename)
        chrom = 'chr16'
        real_donors = [2138327, 30020798]   # chain +, chain -
        real_acceptors = [2138227, 30020879]   # also present in SJ.out.tab
        star_donors = []
        star_acceptors = []
        for real_donor in real_donors:
            self.assertIn(real_donor, donors[chrom])
        for real_acceptor in real_acceptors:
            self.assertIn(real_acceptor, acceptors[chrom])


    def test_get_interval_length(self):
        intervals = [
            interval([100.0, 119.0]),
            interval([100.0, 119.0], [152.0, 171.0]),
            interval([100.0, 119.0], [152.0, 171.0], [1172.0, 1247.0]),
        ]
        exp_lengths = [
            20,
            40,
            116,
        ]

        for i, one_interval in enumerate(intervals):
            res_length = ptes.get_interval_length(one_interval)
            self.assertEqual(exp_lengths[i], res_length)


    def test_dict_to_interval(self):
        exp_intervals = [
            interval([100.0, 119.0], [120.0, 151.0], [152.0, 171.0], [1172.0, 1247.0]),
            interval([100.0, 119.0], [120.0, 151.0], [152.0, 171.0], [1248.0, 1260.0]),
            interval([100.0, 119.0], [120.0, 151.0], [152.0, 171.0], [1248.0, 1260.0]),
        ]
        for i, read_dict in enumerate(self.input_read_dicts):
            res_interval = ptes.dict_to_interval(read_dict=read_dict)
            self.assertEqual(exp_intervals[i], res_interval)


    def test_mate_intersection(self):
        res_list = [
            ['inside', 'outside', 'outside'],
            ['inside', 'inside', 'outside'],
            ['inside', 'inside', 'outside'],
                    ]
        mates1 = [
            interval([100.0, 119.0], [152.0, 171.0]),
            interval([100.0, 119.0], [152.0, 171.0], [1248.0, 1260.0]),
            interval([100.0, 119.0], [120.0, 151.0], [152.0, 171.0], [1248.0, 1260.0])
        ]
        mates2 = [
            interval([120.0, 151.0]),
            interval([1172.0, 1247.0]),
            interval([2000.0, 2047.0]),
        ]
        for i, mate1 in enumerate(mates1):
            for j, mate2 in enumerate(mates2):
                self.assertEqual(ptes.mate_intersection(mate1, mate2), res_list[i][j])


    def test_return_mates(self):
        args = ['cigar1', 'coord1', 'cigar2', 'coord2', 'chain']
        cases = [
            ['100M200p50M20S', 4000, '50S20M', 1000, '+'],
            ['20S50M200p100M', 1000, '20M50S', 2000, '-'],
            ['100M2000p20M50S', 1000, '20S50M', 2000, '+'],
            ['20S50M2000p100M', 1000, '20M50S', 2000, '-'],
            ['20M50S', 2000, '20S50M2000p100M', 1000, '+'],
            ['20S50M', 2000, '100M1900p20M50S', 1000, '-'],
            ['20S50M3000S', 2000, '100M1900p20M50S', 1000, '-'],
            ['20S50M', 2000, '100M910N50M940p20M50S', 1000, '-'],
        ]

        exp_mates = [
            (interval([1000.0, 4349.0]), interval([4000.0, 4099.0])),
            (interval([1000.0, 2019.0]), interval([1250.0, 1349.0])),
            (interval([2000.0, 3119.0]), interval([1000.0, 1099.0])),
            (interval([1000.0, 2019.0]), interval([3050.0, 3149.0])),
            (interval([1000.0, 2019.0]), interval([3050.0, 3149.0])),
            (interval([2000.0, 3019.0]), interval([1000.0, 1099.0])),
            (interval([2000.0, 3019.0]), interval([1000.0, 1099.0])),
            (interval([2000.0, 3019.0]), interval([1000.0, 2059.0])),
            ]
        for i, case in enumerate(cases):
            mate1, mate2 = ptes.return_mates(**dict(zip(args, case)))
            self.assertEqual((mate1,mate2), exp_mates[i])
            print ptes.mate_intersection(mate1, mate2)


    def test_order_interval_list(self):
        interval_lists = [
            [interval([1000.0, 4349.0]), interval([4000.0, 4099.0])],
            [interval([1000.0, 4349.0]), interval([4000.0, 4099.0]), interval([3050.0, 3149.0])],
        ]


if __name__ == "__main__":
    unittest.main()