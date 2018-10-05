import os
import shutil
import unittest

from ptes.lib import general
from ptes.constants import TEST_DIR, PTES_logger


PACKAGE_DIR = 'lib'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestGeneral(unittest.TestCase):

    def test_sign(self):
        exp_sign = -1
        res_sign = general.sign(-10)
        self.assertEqual(exp_sign, res_sign)


if __name__ == "__main__":
    unittest.main()