# -*- coding: utf-8 -*-

import logging

from rtspm.version import __version__

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = [
    '__version__'
]