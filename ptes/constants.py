import os
import logging

from .lib.general import setup_logging


CONSTANTS_PATH = os.path.dirname(os.path.realpath(__file__))
CONF_DIR_NAME = 'configs'
CONF_DIR_PATH = os.path.join(CONSTANTS_PATH, CONF_DIR_NAME)
DEFAULT_CONFIG = os.path.join(CONF_DIR_PATH, "default_config.ini")
LOG_CONFIG_NAME = 'log_conf.yaml'
LOGGER_NAME = 'PTES_logger'

setup_logging(DEFAULT_CONFIG)
PTES_logger = logging.getLogger(LOGGER_NAME)