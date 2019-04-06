import os
import shutil
import unittest
import json
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

from ptes.star import junctions_to_bed
from ptes.constants import TEST_DIR, PTES_logger
from ptes import ptes


PACKAGE_DIR = 'ptes'
INPUT_DIR = os.path.join(TEST_DIR, 'test_data', PACKAGE_DIR)
OUTPUT_DIR = os.path.join(TEST_DIR, 'test_results', PACKAGE_DIR)

try:
    os.makedirs(OUTPUT_DIR)
except OSError:
    pass


class TestJunc(unittest.TestCase):
    def test_main(self):
        junctions_to_bed.main()