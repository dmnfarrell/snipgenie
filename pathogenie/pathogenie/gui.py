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

class App(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, filenames=[], project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("pathogenie")

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
        #mainlayout.addWidget(left)
        self.opts = AppOptions(parent=self.m)
        dialog = self.opts.showDialog(left, wrap=1, section_wrap=1)
        left.setFixedWidth(250)

        center = QWidget(self.m)#, orientation=QtCore.Qt.Vertical)
        #mainlayout.addWidget(center)
        l = QVBoxLayout(center)
        self.fastq_table = tables.FilesTable(center, app=self, dataframe=pd.DataFrame())
        l.addWidget(self.fastq_table)

        self.right = right = QWidget(self.m)
        l2 = QVBoxLayout(right)
        #mainlayout.addWidget(right)
        self.right_tabs = QTabWidget(right)
        self.right_tabs.setTabsClosable(True)
        self.right_tabs.tabCloseRequested.connect(self.close_right_tab)
        l2.addWidget(self.right_tabs)
        self.info = QTextEdit(right, readOnly=True)
        font = QFont("Monospace")
        font.setPointSize(9)
        font.setStyleHint(QFont.TypeWriter)
        self.info.setFont(font)
        self.right_tabs.addTab(self.info, 'log')
        self.info.setText("Welcome")
        self.m.setSizes([50,200,150])
        self.m.setStretchFactor(1,0)

        self.statusBar = QStatusBar()
        from . import __version__
        self.projectlabel = QLabel('')
        self.statusBar.addWidget(self.projectlabel, 1)
        self.outdirLabel = QLabel("")
        self.statusBar.addWidget(self.outdirLabel, 1)
        self.progressbar = QProgressBar()
        self.progressbar.setRange(0,1)
        self.statusBar.addWidget(self.progressbar, 2)
        self.setStatusBar(self.statusBar)
        return

    @QtCore.Slot(int)
    def close_right_tab(self, index):
        """Close right tab"""

        name = self.right_tabs.tabText(index)
        if name == 'log':
            return
        self.right_tabs.removeTab(index)
        return

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('&File', self)
        #self.file_menu.addAction('&New', self.newProject,
        #        QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        self.file_menu.addAction('&Add Fastq Files', self.load_fastq_files_dialog,
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

        self.settings_menu = QMenu('&Settings', self)
        self.menuBar().addMenu(self.settings_menu)
        self.settings_menu.addAction('&Set Output Folder', self.set_output_folder)
        #self.settings_menu.addAction('&Add Blast Sequences', self.add_sequences_db)

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&Help', self.online_documentation)
        self.help_menu.addAction('&About', self.about)

    def save_project(self):
        """Save project"""

        if self.proj_file == None:
            self.save_project_dialog()

        filename = self.proj_file
        data={}
        data['inputs'] = self.fastq_table.getDataFrame()
        #data['sheets'] = self.sheets
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
            if not os.path.splitext(filename)[1] == '.proj':
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

        reply=None
        if ask == True:
            reply = QMessageBox.question(self, 'Confirm', "This will clear the current project.\nAre you sure?",
                                        QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.No:
            return False

        self.outputdir = None
        self.sheets = {}
        self.proj_file = None
        self.fastq_table.setDataFrame(pd.DataFrame({'name':[]}))
        #self.tabs.clear()
        self.projectlabel.setText('')
        self.outdirLabel.setText(self.outputdir)
        return

    def load_project(self, filename=None):
        """Load project"""

        self.clear_project()
        data = pickle.load(open(filename,'rb'))
        keys = ['sheets','annotations','outputdir']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        #for s in self.sheets:
            #df = self.sheets[s]['data']
            #kind = self.sheets[s]['kind']
            #self.add_table(s, df, kind)
        ft = self.fastq_table
        ft.setDataFrame(data['inputs'])
        ft.resizeColumns()
        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        self.outdirLabel.setText(self.outputdir)
        return

    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Project', './',
                                        filter="Project Files(*.proj);;All Files(*.*)")
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
        self.load_fastq_table(filenames)
        return

    def load_samples(self):
        """load samples from text file"""

        return

    def load_fastq_table(self, filenames):
        """Append/Load fasta inputs into table"""

        if filenames is None or len(filenames) == 0:
            return

        info = [tools.get_fastq_info(f) for f in filenames]
        new = pd.DataFrame(info)
        df = self.fastq_table.model.df
        if len(df)>0:
            new = pd.concat([df,new],sort=False).reset_index(drop=True)
        self.fastq_table.setDataFrame(new)
        self.fastq_table.resizeColumns()
        return

    def load_fastq_files_dialog(self):
        """Load fasta files"""

        options = QFileDialog.Options()
        filenames, _ = QFileDialog.getOpenFileNames(self, 'Open File', './',
                                                    filter="Fastq Files(*.fq *.fastq *.fastq.gz);;All Files(*.*)")
        if not filenames:
            return
        self.load_fastq_table(filenames)
        return

    def set_output_folder(self):
        """Set the output folder"""

        selected_directory = QFileDialog.getExistingDirectory()
        if selected_directory:
            self.outputdir = selected_directory
        #check it's empty?
        self.outdirLabel.setText(self.outputdir)
        return

    def check_output_folder(self):
        """check if we have an output dir"""

        if self.outputdir == None:
            #QMessageBox.warning(self, 'No output folder set',
            #    'You should set an output folder from the Settings menu')
            self.show_info('You should set an output folder from the Settings menu')
            return 0
        return 1

    def quit(self):
        self.close()

    def closeEvent(self, ce):
        self.quit()

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def online_documentation(self,event=None):
        """Open the online documentation"""

        import webbrowser
        link='https://github.com/dmnfarrell/btbgenie'
        webbrowser.open(link,autoraise=1)
        return

    def about(self):

        from . import __version__
        import matplotlib
        import PySide2
        pandasver = pd.__version__
        pythonver = platform.python_version()
        mplver = matplotlib.__version__
        qtver = PySide2.QtCore.__version__
        if self._check_snap == True:
            snap='(snap)'
        else:
            snap=''

        text='pathogenie GUI\n'\
            +'version '+__version__+snap+'\n'\
            +'Copyright (C) Damien Farrell 2020-\n'\
            +'This program is free software; you can redistribute it and/or\n'\
            +'modify it under the terms of the GNU GPL\n'\
            +'as published by the Free Software Foundation; either\n'\
            +'version 3 of the License, or (at your option) any\n'\
            +'later version.\n'\
            +'Using Python v%s, PySide2 v%s\n' %(pythonver, qtver)\
            +'pandas v%s, matplotlib v%s' %(pandasver,mplver)

        msg = QMessageBox.about(self, "About", text)
        return

#https://www.learnpyqt.com/courses/concurrent-execution/multithreading-pyqt-applications-qthreadpool/
class Worker(QtCore.QRunnable):
    """Worker thread for running background tasks."""

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.kwargs['progress_callback'] = self.signals.progress

    @QtCore.Slot()
    def run(self):
        try:
            result = self.fn(
                *self.args, **self.kwargs,
            )
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

class WorkerSignals(QtCore.QObject):
    """
    Defines the signals available from a running worker thread.
    Supported signals are:
    finished
        No data
    error
        `tuple` (exctype, value, traceback.format_exc() )
    result
        `object` data returned from processing, anything
    """
    finished = QtCore.Signal()
    error = QtCore.Signal(tuple)
    result = QtCore.Signal(object)
    progress = QtCore.Signal(str)

class AppOptions(widgets.BaseOptions):
    """Class to provide a dialog for global plot options"""

    def __init__(self, parent=None):
        """Setup variables"""
        self.parent = parent
        self.kwds = {}
        aligners = ['bwa','bowtie','bowtie2']
        cpus = [str(i) for i in range(1,os.cpu_count()+1)]
        self.groups = {'general':['threads','overwrite'],
                        'aligners':['aligner'],
                       'blast':['db','identity','coverage'],
                       }
        self.opts = {'threads':{'type':'combobox','default':4,'items':cpus},
                    'overwrite':{'type':'checkbox','default':True},
                    'aligner':{'type':'combobox','default':'bwa',
                    'items':aligners,'label':'aligner'},
                    'db':{'type':'combobox','default':'card',
                    'items':[],'label':'database'},
                    'identity':{'type':'entry','default':90},
                    'coverage':{'type':'entry','default':50}
                    }
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='pathogenie gui tool')
    parser.add_argument("-f", "--fasta", dest="filenames",default=[],
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-p", "--proj", dest="project",default=None,
                        help="load .proj project file", metavar="FILE")
    args = vars(parser.parse_args())

    app = QApplication(sys.argv)
    aw = App(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
