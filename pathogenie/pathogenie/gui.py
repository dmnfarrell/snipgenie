#!/usr/bin/env python

"""
    pathogenie GUI.
    Created Jan 2020
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

from __future__ import absolute_import, print_function
import sys,os,traceback,subprocess
import glob,platform,shutil
import pickle
import threading,time
from PySide2 import QtCore
from PySide2.QtWidgets import *
from PySide2.QtGui import *

import pandas as pd
import numpy as np
from Bio import SeqIO
from . import tools, app, widgets, tables

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
logoimg = os.path.join(module_path, 'logo.png')

class myApp(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, filenames=[], project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("pygenefinder")

        self.setWindowIcon(QIcon(logoimg))
        self.create_menu()
        self.main = QSplitter(self)
        #screen_resolution = QDesktopWidget().screenGeometry()
        screen = QGuiApplication.screens()[0]
        screen_resolution  = screen.geometry()
        if screen_resolution.height() > 1080:
            fac=0.8
            width, height = screen_resolution.width()*fac, screen_resolution.height()*fac
            self.setGeometry(QtCore.QRect(150, 150, width, height))
        else:
            self.showMaximized()

        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.setup_gui()
        self.clear_project()

        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
        if project != None:
            self.load_project(project)
        return
        
    def setup_gui(self):
        """Add all GUI elements"""

        self.m = QSplitter(self.main)
        #mainlayout = QHBoxLayout(self.m)
        left = QWidget(self.m)
        center = QWidget(self.m)#, orientation=QtCore.Qt.Vertical)
        #mainlayout.addWidget(center)
        l = QVBoxLayout(center)
        self.fasta_table = tables.FilesTable(center, app=self, dataframe=pd.DataFrame())
        l.addWidget(self.fasta_table)
        self.fasta_table.setColumnWidth(0,200)
        self.fasta_table.setColumnWidth(1,400)

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('&File', self)
        #self.file_menu.addAction('&New', self.newProject,
        #        QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        self.file_menu.addAction('&Add Fasta Files', self.load_fasta_files_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.file_menu.addAction('&Load Test Files', self.load_test,
                QtCore.Qt.CTRL + QtCore.Qt.Key_T)
        self.menuBar().addSeparator()
        self.file_menu.addAction('&New Project', self.new_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        self.file_menu.addAction('&Open Project', self.load_project_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        self.file_menu.addAction('&Save Project', self.save_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction('&Save Project As', self.save_project_dialog)
        self.file_menu.addAction('&Quit', self.quit,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)
        
    def save_project(self):
        """Save project"""

        if self.proj_file == None:
            self.save_project_dialog()

        filename = self.proj_file
        data={}
        data['inputs'] = self.fasta_table.getDataFrame()
        data['sheets'] = self.sheets
        data['annotations'] = self.annotations
        data['outputdir'] = self.outputdir
        #data['meta'] = self.saveMeta(table)
        self.projectlabel.setText(filename)
        pickle.dump(data, open(filename,'wb'))
        return

    def save_project_dialog(self):
        """Save as project"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Project files (*.proj);;All files (*.*)",
                                                  options=options)
        if filename:
            if not os.path.splitext(filename)[1] == '.pygf':
                filename += '.pygf'
            self.proj_file = filename
            self.save_project()
        return

    def new_project(self):
        """New project"""

        self.clear_project(ask=True)
        #self.set_output_folder()
        return

    def clear_project(self, ask=False):
        """Clear all loaded inputs and results"""
        
        return
        
    def load_project(self, filename=None):
        """Load project"""

        self.clear_project()
        data = pickle.load(open(filename,'rb'))
        
    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Project', './',
                                        filter="Project Files(*.pygf);;All Files(*.*)")
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def load_test(self):
        """Load test_files"""

        reply = self.clear_project(ask=True)
        if reply == False:
            return
        filenames = glob.glob(os.path.join(app.datadir, '*.fa'))
        self.load_fasta_table(filenames)
        return
        
