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
from . import tools, aligners, app, widgets, tables, plotting

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
        app.copy_ref_genomes()
        self.update_ref_genomes()
        self.clear_project()

        if platform.system() == 'Windows':
            app.fetch_binaries()
        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
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

        self.tabs = QTabWidget(center)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        l.addWidget(self.tabs)

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
    def close_tab(self, index):
        """Close current tab"""

        #index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        self.tabs.removeTab(index)
        #del self.sheets[name]
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

        self.analysis_menu = QMenu('&Analysis', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.analysis_menu)
        #self.analysis_menu.addAction('&Run Blast',
        #    lambda: self.run_threaded_process(self.run_gene_finder, self.find_genes_completed))
        self.analysis_menu.addAction('&Trim Reads',
            lambda: self.run_threaded_process(self.run_trimming, self.processing_completed))
        self.analysis_menu.addAction('&Align Reads',
            lambda: self.run_threaded_process(self.align_files, self.processing_completed))
        self.analysis_menu.addAction('&Call Variants',
            lambda: self.run_threaded_process(self.variant_calling, self.processing_completed))
        self.analysis_menu.addAction('&Show Annotation', self.show_ref_annotation)

        self.settings_menu = QMenu('&Settings', self)
        self.menuBar().addMenu(self.settings_menu)
        self.settings_menu.addAction('&Set Output Folder', self.set_output_folder)
        self.settings_menu.addAction('&Add Reference Sequence', self.add_reference)

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
        data['results'] = self.results
        data['inputs'] = self.fastq_table.getDataFrame()
        #data['sheets'] = self.sheets
        data['outputdir'] = self.outputdir

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
                filename += '.proj'
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
        self.results = {}
        self.proj_file = None
        self.fastq_table.setDataFrame(pd.DataFrame({'name':[]}))
        #self.tabs.clear()
        self.ref_genome = app.ref_genome
        self.projectlabel.setText('')
        self.outdirLabel.setText(self.outputdir)
        return

    def load_project(self, filename=None):
        """Load project"""

        self.clear_project()
        data = pickle.load(open(filename,'rb'))
        keys = ['sheets','outputdir','results']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        ft = self.fastq_table
        ft.setDataFrame(data['inputs'])
        ft.resizeColumns()
        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        self.outdirLabel.setText(self.outputdir)
        print (self.results)
        if 'vcf_file' in self.results:
            self.show_variants()
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

    def load_fastq_table(self, filenames):
        """Append/Load fasta inputs into table"""

        if filenames is None or len(filenames) == 0:
            return

        df = self.fastq_table.model.df

        new = app.get_sample_names(filenames)
        new['read_length'] = new.filename.apply(tools.get_fastq_info)
        #print (new)

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

    def add_reference(self):
        """Add a reference genome sequence from fasta"""

        path = app.sequencedir
        if not os.path.exists(path):
            os.makedirs(path)
        filename, _ = QFileDialog.getOpenFileName(self, 'Open Fasta', './',
                                    filter="Fasta Files(*.fa *.fna *.fasta);;All Files(*.*)")
        if not filename:
            return
        #copy file
        dest = os.path.join(path, os.path.basename(filename))
        shutil.copy(filename, dest)
        #update widget
        self.update_ref_genomes()
        return

    def update_ref_genomes(self):

        path = app.sequencedir
        files = glob.glob(os.path.join(path,'*.f*a'))
        labels = [os.path.basename(i) for i in files]
        self.opts.widgets['refgenome'].addItems(labels)
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

    def run_trimming(self, progress_callback):
        """Run quality and adapter trimming"""

        retval = self.check_output_folder()
        if retval == 0:
            return
        progress_callback.emit('Running trimming')
        self.opts.applyOptions()
        kwds = self.opts.kwds
        self.running = True
        overwrite = kwds['overwrite']
        threshold = kwds['quality']
        df = self.fastq_table.model.df
        path = os.path.join(self.outputdir, 'trimmed')
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        st=time.time()
        for i,row in df.iterrows():
            outfile = os.path.join(path, os.path.basename(row.filename))
            progress_callback.emit(outfile)
            if not os.path.exists(outfile) or overwrite == True:
                tools.trim_reads(row.filename, outfile)
            df.loc[i,'trimmed'] = outfile
            self.fastq_table.refresh()
        t = round(time.time()-st,1)
        progress_callback.emit('took %s seconds' %str(t))
        return

    def align_files(self, progress_callback):
        """Run gene annotation for input files.
        progress_callback: signal for indicating progress in gui
        """

        retval = self.check_output_folder()
        if retval == 0:
            return
        self.running = True
        self.opts.applyOptions()
        kwds = self.opts.kwds
        overwrite = kwds['overwrite']
        df = self.fastq_table.model.df

        #rows = self.fastq_table.getSelectedRows()
        #df = df.loc[rows]
        msg = 'Aligning reads..\nThis may take some time.'
        progress_callback.emit(msg)
        ref = app.ref_genome
        aligners.build_bwa_index(ref)

        progress_callback.emit('Using reference genome: %s' %ref)
        path = os.path.join(self.outputdir, 'mapped')
        if not os.path.exists(path):
            os.makedirs(path)
        app.align_reads(df, idx=ref, outdir=path, overwrite=overwrite, threads=kwds['threads'],
                            callback=progress_callback.emit)
        return

    def variant_calling(self, progress_callback):
        """Run variant calling for available bam files."""

        retval = self.check_output_folder()
        if retval == 0:
            return
        self.running = True
        df = self.fastq_table.model.df

        #use trimmed files if present in table
        bam_files = list(df.bam_file.unique())
        print (bam_files)
        path = self.outputdir
        self.results['vcf_file'] = app.variant_calling(bam_files, app.ref_genome,
                                    path, callback=progress_callback.emit)
        self.show_variants()
        return

    def show_variants(self):

        vcf_file = self.results['vcf_file']
        vdf = tools.vcf_to_dataframe(vcf_file)
        table = tables.DefaultTable(self.tabs, app=self, dataframe=vdf)
        i = self.tabs.addTab(table, 'variants')
        self.tabs.setCurrentIndex(i)
        return

    def processing_completed(self):
        """Alignment/calling completed"""

        self.info.append("finished")
        self.progressbar.setRange(0,1)
        df = self.fastq_table.getDataFrame()
        self.fastq_table.refresh()
        self.running = False
        return

    def run(self):
        """Run all steps"""

        self.run_trimming()
        self.align_files()
        self.variant_calling()
        return

    def run_threaded_process(self, process, on_complete):
        """Execute a function in the background with a worker"""

        worker = Worker(fn=process)
        self.threadpool.start(worker)
        worker.signals.finished.connect(on_complete)
        worker.signals.progress.connect(self.progress_fn)
        self.progressbar.setRange(0,0)
        return

    def progress_fn(self, msg):

        print (msg)
        self.info.append(msg)
        self.info.verticalScrollBar().setValue(1)
        return

    def show_ref_annotation(self):

        gff_file = app.ref_gff
        feats = tools.gff_to_features(gff_file)
        df = tools.features_to_dataframe(feats)
        t = tables.DataFrameTable(self.tabs, dataframe=df)
        i = self.tabs.addTab(t, 'ref_annotation')
        self.tabs.setCurrentIndex(i)
        return

    def quality_summary(self, row):
        """Summary of features"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = 'qual:'+data['name']
        if name in self.get_tab_names():
            return
        w = widgets.PlotViewer(self)
        import pylab as plt
        fig,ax = plt.subplots(2,1, figsize=(7,5), dpi=65, facecolor=(1,1,1), edgecolor=(0,0,0))
        axs=ax.flat
        tools.plot_fastq_qualities(data.filename, ax=axs[0])
        tools.plot_fastq_gc_content(data.filename, ax=axs[1])
        plt.tight_layout()
        w.show_figure(fig)
        i = self.tabs.addTab(w, name )
        self.tabs.setCurrentIndex(i)
        return

    def show_bam_viewer(self, row):
        """Show simple alignment view for a bam file"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = 'aln:'+data['name']
        if name in self.get_tab_names():
            return
        w = widgets.BamViewer(self)
        w.load_data(data.bam_file, app.ref_genome)
        w.redraw(start=1, end=2000)
        i = self.tabs.addTab(w, name )
        self.tabs.setCurrentIndex(i)
        return

    def show_info(self, msg):

        self.info.append(msg)
        self.info.verticalScrollBar().setValue(
            self.info.verticalScrollBar().maximum())
        return

    def get_tab_names(self):
        return {self.tabs.tabText(index):index for index in range(self.tabs.count())}

    def quit(self):
        self.close()
        return

    def closeEvent(self, event=None):

        if self.proj_file != None and event != None:
            reply = QMessageBox.question(self, 'Confirm', "Save the current project?",
                                            QMessageBox.Cancel | QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Cancel:
                event.ignore()
                return
            elif reply == QMessageBox.Yes:
                self.save_project()
        #self.close()
        event.accept()

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
            +'This program is free software; you can redistribute it and/or '\
            +'modify it under the terms of the GNU GPL '\
            +'as published by the Free Software Foundation; either '\
            +'version 3 of the License, or (at your option) any '\
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
        genomes = []
        aligners = ['bwa','bowtie','bowtie2']
        cpus = [str(i) for i in range(1,os.cpu_count()+1)]
        self.groups = {'general':['threads','overwrite'],
                        'reference':['refgenome','annotation'],
                        'trimming':['quality'],
                        'aligners':['aligner'],
                       'blast':['db','identity','coverage'],
                       }
        self.opts = {'threads':{'type':'combobox','default':4,'items':cpus},
                    'overwrite':{'type':'checkbox','default':False},
                    'refgenome':{'type':'combobox','default':'',
                    'items':genomes,'label':'genome'},
                    'annotation':{'type':'combobox','default':'',
                    'items':[],'label':'annotation'},
                    'aligner':{'type':'combobox','default':'bwa',
                    'items':aligners,'label':'aligner'},
                    'db':{'type':'combobox','default':'card',
                    'items':[],'label':'database'},
                    'identity':{'type':'entry','default':90},
                    'coverage':{'type':'entry','default':50},
                    'quality':{'type':'spinbox','default':30}
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
