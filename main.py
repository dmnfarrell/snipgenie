#!/usr/bin/env python

"""
    snipgenie GUI app lauancher.
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

try:
    from PySide2.QtWidgets import *
    from PySide2 import QtCore
except:
    from PyQt5.QtWidgets import *
    from PyQt5 import QtCore
from snipgenie import gui

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snipgenie tool')
    parser.add_argument("-f", "--fasta", dest="filename",
                        help="input fasta file", metavar="FILE")
    args = vars(parser.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = gui.App()
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
