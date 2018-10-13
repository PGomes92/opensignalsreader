# -*- coding: utf-8 -*-
# Compatibility
from __future__ import absolute_import

# allow lazy import
from .opensignalsreader import OpenSignalsReader

# get version
from .__version__ import __version__

# Package info
__author__ = "Pedro Gomes"
__email__ = "pgomes92@gmail.com"
__maintainer__ = "Pedro Gomes"
__status__ = "Development"
__license__ = "MIT"
name = "opensignalsreader"
description = "Python package to read OpenSignals (r)evolution files."
long_description = "Python package to read OpenSignals (r)evolution files and automatic sensor data conversion " \
				   "using BITalino (r)evolution transfer functions."