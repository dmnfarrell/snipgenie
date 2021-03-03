#!/usr/bin/env python
"""
    Implements tests for snpgenie
    Created Sep 2020
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys, os, tempfile
from . import app, tools, aligners, trees
import unittest
tempdir = tempfile.gettempdir()
module_path = os.path.dirname(os.path.abspath(__file__))
testdatadir = os.path.join(module_path, 'testing')

class snpgenieTests(unittest.TestCase):
    """snpegenie tests"""
    def setUp(self):
        return

    def test_workflow(self):
        """Workflow run test"""

        out = os.path.join(tempdir, 'snpgenie_tests')
        args = {'threads':8, 'outdir': out, 'input': testdatadir,
                'aligner':'bowtie',
                'reference': None, 'overwrite':False}
        W = app.WorkFlow(**args)
        st=W.setup()
        if st == True:
            W.run()
        return


if __name__ == '__main__':
    unittest.main()
