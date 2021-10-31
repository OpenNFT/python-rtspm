# -*- coding: utf-8 -*-

from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version('rtspm')
except PackageNotFoundError:  # pragma: no cover
    __version__ = '0.0.0.dev'
