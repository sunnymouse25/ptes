import os
import shutil
import unittest


from ptes.panhandles import random_shuffle
from ptes.constants import TEST_DIR, PTES_logger


PACKAGE_DIR = 'ptes'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestPanhandles(unittest.TestCase):
    def test_choose_close(self):
        sorted_list = range(100)
        res_dict = random_shuffle.choose_close(sorted_list=sorted_list,
                                               threshold=5)
        for key, value in res_dict.items():
            print key, value
