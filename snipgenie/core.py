#!/usr/bin/env python

"""
    snipgenie methods for cmd line tool.
    Created Nov 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warroanty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,subprocess,glob,re
import time, datetime
import platform
import pandas as pd
import numpy as np

if 'Windows' in platform.platform():
    defaultfont = 'Arial'
else:
    defaultfont = 'Monospace'
defaults = {
            'FONT' :defaultfont,
            'FONTSIZE' : 10,
            'TIMEFORMAT' :'%m/%d/%Y',
            'PLOTSTYLE' : 'bmh',
            'ICONSIZE' : 20,
            'DPI' : 100,
            'BGCOLOR' : '#F4F4F3',
            'THEME': 'Fusion'
}
#populate current class variable
for k in defaults:
    vars()[k] = defaults[k]
