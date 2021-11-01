#!/usr/bin/env python

"""
    snipgenie GUI.
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

import sys,os,traceback,subprocess
import glob,platform,shutil
import pickle
import threading,time
from .qt import *
import pandas as pd
import numpy as np
import pylab as plt
from Bio import SeqIO
from . import tools, aligners, app, widgets, tables, plotting, trees

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
logoimg = os.path.join(module_path, 'logo.png')
stylepath = os.path.join(module_path, 'styles')
iconpath = os.path.join(module_path, 'icons')

class App(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, filenames=[], project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("snipgenie")

        self.setWindowIcon(QIcon(logoimg))
        self.create_menu()
        self.main = QSplitter(self)
        screen = QGuiApplication.screens()[0]
        screen_resolution  = screen.geometry()
        if screen_resolution.height() > 1280:
            fac=0.8
            width, height = screen_resolution.width()*fac, screen_resolution.height()*fac
            self.setGeometry(QtCore.QRect(150, 150, width, height))
        else:
            self.showMaximized()

        self.recent_files = ['']
        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.create_tool_bar()
        self.setup_gui()
        self.load_settings()
        self.show_recent_files()

        self.new_project()
        self.running = False

        if platform.system() == 'Windows':
            app.fetch_binaries()
        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
        return

    def load_settings(self):
        """Load GUI settings"""

        s = self.settings = QtCore.QSettings('snipgenie','default')
        try:
            self.resize(s.value('window_size'))
            self.move(s.value('window_position'))
            self.setStyle(s.value('style'))
            #self.FONT = s.value("font")
            #self.FONTSIZE = int(s.value("fontsize"))
            r = s.value("recent_files")
            if r != '':
                rct = r.split(',')
                self.recent_files = [f for f in rct if os.path.exists(f)]

        except Exception as e:
            print (e)
        return

    def save_settings(self):
        """Save GUI settings"""

        self.settings.setValue('window_size', self.size())
        self.settings.setValue('window_position', self.pos())
        #self.settings.setValue('style', self.style)
        #self.settings.setValue('font', self.FONT)
        #self.settings.setValue('fontsize', self.FONTSIZE)
        self.settings.setValue('recent_files',','.join(self.recent_files))
        self.settings.sync()
        return

    def set_style(self, style='default'):
        """Change interface style."""

        if style == 'default':
            	self.setStyleSheet("")
        else:
            f = open(os.path.join(stylepath,'%s.qss' %style), 'r')
            self.style_data = f.read()
            f.close()
            self.setStyleSheet(self.style_data)
        self.style = style
        return

    def update_ref_genome(self):
        """Update the ref genome labels"""

        if self.ref_genome == None:
            refname=''
        else:
            refname = os.path.basename(self.ref_genome)
        if self.ref_gb == None:
            gbname = ''
        else:
            gbname = os.path.basename(self.ref_gb)
        self.reflabel.setText(refname)
        self.annotlabel.setText(gbname)
        return

    def update_mask(self):
        if self.mask_file == None:
            self.masklabel.setText('')
        else:
            self.masklabel.setText(os.path.basename(self.mask_file))
        return

    def create_tool_bar(self):
        """Create main toolbar"""

        items = {'New project': {'action': lambda: self.new_project(ask=True),'file':'document-new'},
                 'Open': {'action':self.load_project,'file':'document-open'},
                 'Save': {'action': lambda: self.save_project(),'file':'save'},
                 'Zoom out': {'action':self.zoom_out,'file':'zoom-out'},
                 'Zoom in': {'action':self.zoom_in,'file':'zoom-in'},
                 'Align Reads': {
                    'action': lambda: self.run_threaded_process(self.align_files, self.alignment_completed),
                    'file':'align-reads'},
                 'Call Variants': {
                    'action': lambda: self.run_threaded_process(self.variant_calling, self.processing_completed),
                    'file':'call-variants'},
                 'Make Phylogeny': {
                    'action': lambda: self.run_threaded_process(self.make_phylo_tree, self.phylogeny_completed),
                    'file':'phylogeny'},
                 'Quit': {'action':self.quit,'file':'application-exit'}
                }

        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
        for i in items:
            if 'file' in items[i]:
                iconfile = os.path.join(iconpath,items[i]['file']+'.png')
                icon = QIcon(iconfile)
            else:
                icon = QIcon.fromTheme(items[i]['icon'])
            btn = QAction(icon, i, self)
            btn.triggered.connect(items[i]['action'])
            #btn.setCheckable(True)
            toolbar.addAction(btn)
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.m = QSplitter(self.main)
        #left menu
        left = QWidget()
        self.m.addWidget(left)
        l = QVBoxLayout(left)
        lbl = QLabel("Reference Genome:")
        l.addWidget(lbl)
        self.reflabel = QLabel()
        self.reflabel.setStyleSheet('color: blue')
        l.addWidget(self.reflabel)
        lbl = QLabel("Reference Annotation:")
        l.addWidget(lbl)
        self.annotlabel = QLabel()
        self.annotlabel.setStyleSheet('color: blue')
        l.addWidget(self.annotlabel)
        lbl = QLabel("Mask:")
        l.addWidget(lbl)
        self.masklabel = QLabel()
        self.masklabel.setStyleSheet('color: blue')
        l.addWidget(self.masklabel)

        #create option widgets
        self.opts = AppOptions(parent=self.m)
        dialog = self.opts.showDialog(left, wrap=1, section_wrap=1)
        l.addWidget(dialog)

        l.addStretch()
        left.setFixedWidth(250)

        center = QWidget()
        self.m.addWidget(center)
        l = QVBoxLayout(center)
        self.fastq_table = tables.FilesTable(center, app=self, dataframe=pd.DataFrame())
        l.addWidget(self.fastq_table)

        self.tabs = QTabWidget(center)
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        l.addWidget(self.tabs)

        self.right = right = QWidget()
        self.m.addWidget(self.right)
        l2 = QVBoxLayout(right)
        #mainlayout.addWidget(right)

        self.right_tabs = QTabWidget(right)
        self.right_tabs.setTabsClosable(True)
        self.right_tabs.tabCloseRequested.connect(self.close_right_tab)
        l2.addWidget(self.right_tabs)
        self.info = widgets.Editor(right, readOnly=True, fontsize=11)
        self.right_tabs.addTab(self.info, 'log')
        self.info.append("Welcome\n")
        self.m.setSizes([50,200,150])
        self.m.setStretchFactor(1,0)

        self.statusBar = QStatusBar()
        from . import __version__
        self.projectlabel = QLabel('')
        self.projectlabel.setStyleSheet('color: blue')
        self.projectlabel.setAlignment(Qt.AlignLeft)
        self.statusBar.addWidget(self.projectlabel, 1)
        lbl = QLabel("Output folder:")
        self.statusBar.addWidget(lbl,1)
        lbl.setAlignment(Qt.AlignRight)
        self.outdirLabel = QLabel("")
        self.statusBar.addWidget(self.outdirLabel, 1)
        self.outdirLabel.setAlignment(Qt.AlignLeft)
        self.outdirLabel.setStyleSheet('color: blue')

        self.progressbar = QProgressBar()
        self.progressbar.setRange(0,1)
        self.statusBar.addWidget(self.progressbar, 3)
        self.progressbar.setAlignment(Qt.AlignRight)
        self.setStatusBar(self.statusBar)

        #redirect stdout
        self._stdout = StdoutRedirect()
        self._stdout.start()
        self._stdout.printOccur.connect(lambda x : self.info.insert(x))
        return

    @QtCore.Slot(int)
    def close_tab(self, index):
        """Close current tab"""

        #index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        self.tabs.removeTab(index)
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
        self.file_menu.addAction('&Add Folder', self.load_fastq_folder_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.file_menu.addAction('&Add Fastq Files', self.load_fastq_files_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_A)
        #self.file_menu.addAction('&Load Test Files', self.load_test,
        #        QtCore.Qt.CTRL + QtCore.Qt.Key_T)
        self.menuBar().addSeparator()
        icon = QIcon(os.path.join(iconpath,'document-new.png'))
        self.file_menu.addAction(icon, '&New Project', lambda: self.new_project(ask=True),
                QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        icon = QIcon(os.path.join(iconpath,'document-open.png'))
        self.file_menu.addAction(icon, '&Open Project', self.load_project_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        self.recent_files_menu = QMenu("Recent Projects",
            self.file_menu)
        self.file_menu.addAction(self.recent_files_menu.menuAction())
        icon = QIcon(os.path.join(iconpath,'save.png'))
        self.file_menu.addAction(icon, '&Save Project', self.save_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction('&Save Project As', self.save_project_dialog)
        icon = QIcon(os.path.join(iconpath,'application-exit.png'))
        self.file_menu.addAction(icon, '&Quit', self.quit,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.view_menu = QMenu('&View', self)
        self.menuBar().addMenu(self.view_menu)
        icon = QIcon(os.path.join(iconpath,'zoom-in.png'))
        self.view_menu.addAction(icon, 'Zoom In', self.zoom_in,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Equal)
        icon = QIcon(os.path.join(iconpath,'zoom-out.png'))
        self.view_menu.addAction(icon, 'Zoom Out', self.zoom_out,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Minus)

        self.style_menu = QMenu("Styles",  self.view_menu)
        self.style_menu.addAction('Default', self.set_style)
        self.style_menu.addAction('Light', lambda: self.set_style('light'))
        self.style_menu.addAction('Dark', lambda: self.set_style('dark'))
        self.view_menu.addAction(self.style_menu.menuAction())

        self.analysis_menu = QMenu('&Analysis', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.analysis_menu)
        self.analysis_menu.addAction('Trim Reads',
            lambda: self.run_threaded_process(self.run_trimming, self.processing_completed))
        icon = QIcon(os.path.join(iconpath,'align-reads.png'))
        self.analysis_menu.addAction(icon, 'Align Reads',
            lambda: self.run_threaded_process(self.align_files, self.alignment_completed))
        icon = QIcon(os.path.join(iconpath,'call-variants.png'))
        self.analysis_menu.addAction(icon, 'Call Variants',
            lambda: self.run_threaded_process(self.variant_calling, self.processing_completed))
        self.analysis_menu.addAction('Create SNP alignment',
            lambda: self.run_threaded_process(self.snp_alignment, self.snp_align_completed))
        self.analysis_menu.addAction('Make Phylogeny',
            lambda: self.run_threaded_process(self.make_phylo_tree, self.phylogeny_completed))
        #self.analysis_menu.addSeparator()
        #self.analysis_menu.addAction('Run Workflow', self.run)

        self.tools_menu = QMenu('Tools', self)
        self.menuBar().addMenu(self.tools_menu)

        self.tools_menu.addAction('Get Read Stats',
            lambda: self.run_threaded_process(self.add_read_lengths, self.processing_completed))
        self.tools_menu.addAction('Get Mapping Stats',
            lambda: self.run_threaded_process(self.add_mapping_stats, self.processing_completed))
        self.tools_menu.addAction('Fastq Qualities Report', self.fastq_quality_report)
        self.tools_menu.addAction('Show Annotation', self.show_ref_annotation)
        self.tools_menu.addAction('Plot SNP Matrix', self.plot_snp_matrix)
        #self.tools_menu.addAction('Map View', self.show_map)
        self.tools_menu.addAction('Tree Viewer', self.tree_viewer)
        self.tools_menu.addSeparator()
        self.tools_menu.addAction('Check Heterozygosity', self.check_heterozygosity)
        self.tools_menu.addAction('RD Analysis (MTBC)',
            lambda: self.run_threaded_process(self.rd_analysis, self.rd_analysis_completed))
        self.tools_menu.addAction('M.bovis Spoligotyping',
            lambda: self.run_threaded_process(self.spoligotyping, self.spotyping_completed))
        #self.tools_menu.addAction('M.bovis SNP typing',
        #    lambda: self.run_threaded_process(self.snp_typing, self.processing_completed))

        self.settings_menu = QMenu('&Settings', self)
        self.menuBar().addMenu(self.settings_menu)
        self.settings_menu.addAction('Set Output Folder', self.set_output_folder)
        self.settings_menu.addAction('Set Reference Sequence', self.set_reference)
        self.settings_menu.addAction('Set Annnotation (genbank)', self.set_annotation)
        #self.settings_menu.addAction('Set Filters', self.set_filters)
        self.settings_menu.addAction('Add Mask File', self.add_mask)
        self.settings_menu.addAction('Add Sample Meta Data', self.merge_meta_data)
        self.settings_menu.addAction('Clean Up Files', self.clean_up)

        self.presets_menu = QMenu('&Load Preset', self)
        self.menuBar().addMenu(self.presets_menu)
        self.load_presets_menu()

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&Help', self.online_documentation)
        self.help_menu.addAction('&About', self.about)

    def load_presets_menu(self):
        """Add preset genomes to menu"""

        genomes = app.preset_genomes
        for name in genomes:
            seqname = genomes[name]['sequence']
            gbfile = genomes[name]['gb']
            if 'mask' in genomes[name]:
                mask = genomes[name]['mask']
            else:
                mask = None
            def func(seqname,gbfile,mask):
                self.set_reference(seqname)
                self.set_annotation(gbfile)
                self.set_mask(mask)
            receiver = lambda seqname=seqname, gb=gbfile, mask=mask: func(seqname, gb, mask)
            self.presets_menu.addAction('%s' %name, receiver)
        return

    def show_recent_files(self):
        """Populate recent files menu"""

        from functools import partial
        if self.recent_files == None:
            return
        for fname in self.recent_files:
            self.recent_files_menu.addAction(fname, partial(self.load_project, fname))
        self.recent_files_menu.setEnabled(len(self.recent_files))
        return

    def add_recent_file(self, fname):
        """Add file to recent if not present"""

        fname = os.path.abspath(fname)
        if fname and fname not in self.recent_files:
            self.recent_files.insert(0, fname)
            if len(self.recent_files) > 5:
                self.recent_files.pop()
        self.recent_files_menu.setEnabled(len(self.recent_files))
        return

    def save_project(self):
        """Save project"""

        if self.proj_file == None:
            self.save_project_dialog()

        filename = self.proj_file
        data={}
        data['inputs'] = self.fastq_table.getDataFrame()
        keys = ['outputdir','results','ref_genome','ref_gb','mask_file']
        for k in keys:
            if hasattr(self, k):
                data[k] = self.__dict__[k]
        if hasattr(self, 'gisviewer'):
            d=self.gisviewer.saveData()
            data['gisviewer'] = d
        if hasattr(self, 'treeviewer'):
            d=self.treeviewer.saveData()
            data['treeviewer'] = d
        self.opts.applyOptions()
        data['options'] = self.opts.kwds
        self.projectlabel.setText(filename)
        pickle.dump(data, open(filename,'wb'))
        self.add_recent_file(filename)
        return

    def save_project_dialog(self):
        """Save as project"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Project files (*.snipgenie);;All files (*.*)",
                                                  options=options)
        if filename:
            if not os.path.splitext(filename)[1] == '.snipgenie':
                filename += '.snipgenie'
            self.proj_file = filename
            self.save_project()
        return

    def new_project(self, ask=False):
        """Clear all loaded inputs and results"""

        reply=None
        if ask == True:
            reply = QMessageBox.question(self, 'Confirm',
                                "This will clear the current project.\nAre you sure?",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.No:
            return False

        self.outputdir = None
        self.sheets = {}
        self.results = {}
        self.proj_file = None
        self.fastq_table.setDataFrame(pd.DataFrame({'name':[]}))
        self.tabs.clear()
        if hasattr(self, 'treeviewer'):
            self.treeviewer.clear()
        self.ref_genome = None
        self.ref_gb = None
        self.mask_file = None
        self.projectlabel.setText('')
        self.outdirLabel.setText(self.outputdir)
        self.clear_tabs()
        self.update_ref_genome()
        self.update_mask()
        return

    def load_project(self, filename=None):
        """Load project"""

        self.new_project()
        data = pickle.load(open(filename,'rb'))
        keys = ['sheets','outputdir','results','ref_genome','ref_gb','mask_file']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        if 'options' in data:
            self.opts.updateWidgets(data['options'])
        ft = self.fastq_table
        ft.setDataFrame(data['inputs'])
        ft.resizeColumns()
        self.update_ref_genome()
        self.update_mask()
        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        self.outdirLabel.setText(self.outputdir)
        #print (self.results)
        #if 'vcf_file' in self.results:
        #    self.show_variants()
        if 'snp_dist' in self.results:
            self.show_snpdist()
        #load any saved maps
        #if 'gisviewer' in data.keys():
            #self.show_map()
            #self.gisviewer.loadData(data['gisviewer'])
        #load tree view
        if 'treeviewer' in data.keys():
            self.tree_viewer()
            self.treeviewer.loadData(data['treeviewer'])
        self.right_tabs.setCurrentIndex(0)
        self.add_recent_file(filename)
        self.check_fastq_table()
        return

    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Project', './',
                                        filter="Project Files(*.snipgenie);;All Files(*.*)")
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def load_test(self):
        """Load test_files"""

        reply = self.new_project(ask=True)
        if reply == False:
            return
        #filenames = glob.glob(os.path.join(app.datadir, '*.fa'))
        self.load_fastq_table(filenames)
        return

    def import_results_folder(self, path):
        """Import previously made results"""

        df = read.csv(os.path.join(path, 'samples.csv'))
        return

    def check_missing_files(self):
        """Check folders for missing files"""

        folders = ['mapped','trimmed']
        #for f in folders:

        return

    def clear_tabs(self):
        """Clear tabbed panes"""

        #self.right_tabs
        return

    def clean_up(self):
        """Clean up intermediate files"""

        msg = 'This will remove all intermediate files in the output folder. Proceed?'
        reply = QMessageBox.question(self, 'Warning!', msg,
                                        QMessageBox.No | QMessageBox.Yes )
        if reply == QMessageBox.No:
            return

        folders = ['mapped','trimmed']
        for l in folders:
            files = glob.glob(os.path.join(self.outputdir, l, '*'))
            print (files)
            for f in files:
                os.remove(f)
        bcf = os.path.join(self.outputdir, 'raw.bcf')
        if os.path.exists(bcf):
             os.remove(bcf)
        df = self.fastq_table.model.df
        cols = ['bam_file','mapped','total']
        for col in cols:
            if col in df.columns:
                df = df.drop(columns=col)
        self.fastq_table.model.df = df
        self.fastq_table.refresh()
        self.results = {}
        self.treefile = None
        return

    def check_fastq_table(self):
        """Update samples file to reflect table"""

        df = self.fastq_table.model.df
        app.write_samples(df, self.outputdir)
        return

    def load_fastq_table(self, filenames):
        """Append/Load fasta inputs into table"""

        if filenames is None or len(filenames) == 0:
            return
        retval = self.check_output_folder()
        if retval == 0:
            return
        self.opts.applyOptions()
        kwds = self.opts.kwds
        df = self.fastq_table.model.df
        new = app.get_samples(filenames, sep=kwds['labelsep'])
        #pivoted
        new = app.get_pivoted_samples(new)
        new['read_length'] = new.filename1.apply(tools.get_fastq_info)

        if len(df)>0:
            new = pd.concat([df,new],sort=False).reset_index(drop=True)
            new = new.drop_duplicates('filename1')
        self.fastq_table.setDataFrame(new)
        self.fastq_table.resizeColumns()
        app.write_samples(df[['sample']], self.outputdir)
        return

    def load_fastq_folder_dialog(self):
        """Load fastq folder"""

        path = QFileDialog.getExistingDirectory(self, 'Add Input Folder', './')
        if not path:
            return
        filenames = app.get_files_from_paths(path)
        self.load_fastq_table(filenames)
        return

    def load_fastq_files_dialog(self):
        """Load fastq files"""

        filenames, _ = QFileDialog.getOpenFileNames(self, 'Open File', './',
                                                    filter="Fastq Files(*.fq *.fastq *.fq.gz *.fastq.gz);;All Files(*.*)")
        if not filenames:
            return
        self.load_fastq_table(filenames)
        return

    def set_reference(self, filename=None):
        """Reset the reference sequence"""

        msg = "This will change the reference genome. You will need to re-run any previous alignments. Are you sure?"
        reply = QMessageBox.question(self, 'Warning!', msg,
                                        QMessageBox.No | QMessageBox.Yes )
        if reply == QMessageBox.No:
            return
        if filename == None:
            filter="Fasta Files(*.fa *.fna *.fasta)"
            filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="%s;;All Files(*.*)" %filter)
            if not filename:
                return
        self.ref_genome = filename
        self.update_ref_genome()
        return

    def set_annotation(self, filename=None):

        if filename == None:
            filename = self.add_file("Genbank Files(*.gb *.gbk *.gbff)")
        self.ref_gb = filename
        #put annotation in a dataframe
        self.annot = tools.genbank_to_dataframe(self.ref_gb)
        self.update_ref_genome()
        return

    def set_mask(self, filename):

        self.mask_file = filename
        self.update_mask()

    def add_mask(self, filename):
        """Add mask bed file"""

        msg = "This will add a bed file as a mask. See help for format. Continue?"
        reply = QMessageBox.question(self, 'Add bed file', msg,
                                        QMessageBox.No | QMessageBox.Yes )
        if reply == QMessageBox.No:
            return
        filename = self.add_file("Bed files(*.bed)")
        self.mask_file = filename
        self.update_mask()
        return

    def merge_meta_data(self):
        """Add sample meta data by merging with file table"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                    filter="Text Files(*.csv *.txt *.tsv);;All Files(*.*)")
        if not filename:
            return

        meta = pd.read_csv(filename)
        dlg = widgets.MergeDialog(self, self.fastq_table, meta)
        dlg.exec_()
        if not dlg.accepted:
            return
        return

    def add_file(self, filter="Fasta Files(*.fa *.fna *.fasta)", path=None):
        """Add a file to the config folders"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                    filter="%s;;All Files(*.*)" %filter)
        if not filename:
            return
        if path != None:
            #copy file
            if not os.path.exists(path):
                os.makedirs(path)
            dest = os.path.join(path, os.path.basename(filename))
            shutil.copy(filename, dest)
            self.show_info('added %s' %filename)
        return filename

    def set_output_folder(self):
        """Set the output folder"""

        selected_directory = QFileDialog.getExistingDirectory()
        if selected_directory:
            self.outputdir = selected_directory
        #check if folder already got some results
        results_file = os.path.join(self.outputdir, 'samples.csv')
        if os.path.exists(results_file):
            msg = "This folder appears to have results already. Try to import them?"
            reply = QMessageBox.question(self, 'Confirm', msg,
                                        QMessageBox.Cancel | QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Cancel:
                return
            elif reply == QMessageBox.Yes:
                df = pd.read_csv(results_file)
                self.fastq_table.model.df = df
                self.fastq_table.refresh()
                self.results['vcf_file'] = os.path.join(self.outputdir, 'filtered.vcf.gz')
        self.outdirLabel.setText(self.outputdir)
        return

    def check_output_folder(self):
        """Check if we have an output dir"""

        if self.outputdir == None:
            self.show_info('You should set an output folder from the Settings menu!',
                        color='red')
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
        if self.ref_genome == None:
            print('no reference genome!')
            return
        self.running = True
        self.opts.applyOptions()
        kwds = self.opts.kwds
        overwrite = kwds['overwrite']
        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        new = df.loc[rows]
        #print (new)
        print ('Aligning reads. This may take some time.')
        ref = self.ref_genome
        if kwds['aligner'] == 'bwa':
            aligners.build_bwa_index(self.ref_genome)
        elif kwds['aligner'] == 'subread':
            aligners.build_subread_index(self.ref_genome)

        print('Using reference genome: %s' %ref)
        path = os.path.join(self.outputdir, 'mapped')
        if not os.path.exists(path):
            os.makedirs(path)
        new = app.align_reads(new, idx=ref, outdir=path, overwrite=overwrite,
                        threads=kwds['threads'],
                        aligner=kwds['aligner'],
                        platform=kwds['platform'],
                        callback=progress_callback.emit)
        self.update_table(new)
        #df.to_csv(os.path.join(self.outputdir,'samples.csv'),index=False)
        #summ = app.results_summary(samples)
        #summ.to_csv(os.path.join(self.outputdir,'summary.csv'),index=False)
        #rewrite samples in case check_missing
        #app.write_samples(df, self.outputdir)
        return

    def update_table(self, new):
        """Update table with changed rows"""

        df = self.fastq_table.model.df
        print (new)
        cols = new.columns
        df.loc[df['sample'].isin(new['sample']), cols] = new[cols]
        return

    def variant_calling(self, progress_callback=None):
        """Run variant calling for available bam files."""

        retval = self.check_output_folder()
        if retval == 0:
            return
        self.running = True
        self.opts.applyOptions()
        kwds = self.opts.kwds
        overwrite = kwds['overwrite']
        threads = int(kwds['threads'])
        filters = kwds['filters']
        proximity = kwds['proximity']
        df = self.fastq_table.model.df
        path = self.outputdir

        gff_file = os.path.join(path, self.ref_gb+'.gff')
        tools.gff_bcftools_format(self.ref_gb, gff_file)

        bam_files = list(df.bam_file.dropna().unique())
        print ('Calling with %s of %s samples aligned' %(len(bam_files), len(df)))
        self.results['vcf_file'] = app.variant_calling(bam_files, self.ref_genome, path,
                                    threads=threads, relabel=True,
                                    overwrite=overwrite, filters=filters,
                                    gff_file=gff_file, mask=self.mask_file,
                                    custom_filters=proximity,
                                    callback=progress_callback.emit)
        self.results['nuc_matrix'] = os.path.join(self.outputdir, 'core.txt')
        self.results['csq_matrix'] = os.path.join(self.outputdir, 'csq.matrix')
        return

    def show_variants(self):
        """Show the stored results from variant calling as tables"""

        vcf_file = self.results['vcf_file']
        vdf = tools.vcf_to_dataframe(vcf_file)
        table = tables.DefaultTable(self.tabs, app=self, dataframe=vdf)
        i = self.tabs.addTab(table, 'variants')
        if 'nuc_matrix' in self.results:
            nucmat = pd.read_csv(self.results['nuc_matrix'],sep=' ')
            table = tables.DefaultTable(self.tabs, app=self, dataframe=nucmat)
            i = self.tabs.addTab(table, 'snp_matrix')
            self.tabs.setCurrentIndex(i)
        if 'csq_matrix' in self.results:
            csqmat = pd.read_csv(self.results['csq_matrix'])
            table = tables.DefaultTable(self.tabs, app=self, dataframe=csqmat)
            i = self.tabs.addTab(table, 'csq_matrix')
            self.tabs.setCurrentIndex(i)
        return

    def show_snpdist(self):
        """Show SNP distance matrix"""

        filename = self.results['snp_dist']
        if not os.path.exists(filename):
            return
        mat = pd.read_csv(filename)
        table = tables.DefaultTable(self.tabs, app=self, dataframe=mat)
        i = self.tabs.addTab(table, 'snp_dist')
        return

    def snp_alignment(self, progress_callback):
        """Make snp matrix from variant positions"""

        self.opts.applyOptions()
        kwds = self.opts.kwds
        vcf_file = self.results['vcf_file']
        print('Making SNP alignment')
        result, smat = tools.fasta_alignment_from_vcf(vcf_file,
                                                callback=progress_callback.emit)
        #print (result)
        outfasta = os.path.join(self.outputdir, 'core.fa')
        self.results['snp_file'] = outfasta
        SeqIO.write(result, outfasta, 'fasta')
        self.results['snp_dist'] = os.path.join(self.outputdir, 'snpdist.csv')
        from Bio import AlignIO
        aln = AlignIO.read(outfasta, 'fasta')
        snp_dist = tools.snp_dist_matrix(aln)
        snp_dist.to_csv(self.results['snp_dist'], sep=',')
        return

    def snp_align_completed(self):

        self.processing_completed()
        self.show_snpdist()
        return

    def make_phylo_tree(self, progress_callback=None, method='raxml'):
        """Make phylogenetic tree"""

        corefasta = os.path.join(self.outputdir, 'core.fa')
        bootstraps = 100
        if method == 'raxml':
            outfile = os.path.join(self.outputdir,'RAxML_bipartitions.variants')
            treefile = trees.run_RAXML(corefasta, bootstraps=bootstraps, outpath=self.outputdir)
        elif method == 'fasttree':
            outfile = os.path.join(self.outputdir,'fasttree.newick')
            treefile = trees.run_fasttree(corefasta, bootstraps=bootstraps, outpath=self.outputdir)

        self.treefile = outfile
        return

    def phylogeny_completed(self):

        self.processing_completed()
        self.show_tree()

    def show_tree(self):
        """Show current tree"""

        self.tree_viewer()
        filename = os.path.join(self.outputdir,'RAxML_bipartitions.variants')
        self.treeviewer.load_tree(filename)
        self.treeviewer.update()
        return

    def tree_viewer(self):
        """Show tree viewer"""

        from . import phylo
        if not hasattr(self, 'treeviewer'):
            self.treeviewer = phylo.TreeViewer(self)
        if not 'phylogeny' in self.get_tabs():
            idx = self.right_tabs.addTab(self.treeviewer, 'phylogeny')
            self.right_tabs.setCurrentIndex(idx)
        return

    def show_map(self):

        from . import gis
        if not hasattr(self, 'gisviewer'):
            self.gisviewer = gis.GISViewer()
        if not 'map' in self.get_tabs():
            idx = self.right_tabs.addTab(self.gisviewer, 'map')
            self.right_tabs.setCurrentIndex(idx)
        return

    def get_tabs(self):

        n=[]
        for i in range(self.right_tabs.count()):
            n.append(self.right_tabs.tabText(i))
        return n

    def processing_completed(self):
        """Generic process completed"""

        self.progressbar.setRange(0,1)
        self.running = False
        self.fastq_table.refresh()
        print ('finished')
        return

    def alignment_completed(self):
        """Alignment/calling completed"""

        print("finished")
        self.progressbar.setRange(0,1)
        df = self.fastq_table.getDataFrame()
        self.fastq_table.refresh()
        self.running = False
        return

    def run(self):
        """Run all steps"""

        #self.run_threaded_process(self.run_trimming, self.processing_completed)
        return

    def run_threaded_process(self, process, on_complete):
        """Execute a function in the background with a worker"""

        if self.running == True:
            return
        worker = Worker(fn=process)
        self.threadpool.start(worker)
        worker.signals.finished.connect(on_complete)
        worker.signals.progress.connect(self.progress_fn)
        self.progressbar.setRange(0,0)
        return

    def progress_fn(self, msg):

        self.info.append(msg)
        self.info.verticalScrollBar().setValue(1)
        return

    def show_ref_annotation(self):
        """Show annotation in table"""

        gb_file = self.ref_gb
        df = tools.genbank_to_dataframe(gb_file)
        t = tables.DataFrameTable(self.tabs, dataframe=df)
        i = self.tabs.addTab(t, 'ref_annotation')
        self.tabs.setCurrentIndex(i)
        return

    def sample_details(self, row):

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        pd.set_option('display.max_colwidth', None)
        print (data)
        print ()
        return

    def quality_summary(self, row):
        """Summary of a single sample, both files"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        if 'filename2' in df.columns:
            colnames = zip([data.name1, data.name2], [data.filename1, data.filename2])
        else:
            colnames = [[data.name1, data.filename1]]
        for name,fname in colnames:
            label = 'qual:'+name
            if label in self.get_tab_names():
                return
            w = widgets.PlotViewer(self)
            fig,ax = plt.subplots(2,1, figsize=(7,5), dpi=65, facecolor=(1,1,1), edgecolor=(0,0,0))
            axs=ax.flat
            if not os.path.exists(fname):
                self.show_info('This file is missing.')
                return
            tools.plot_fastq_qualities(fname, ax=axs[0])
            tools.plot_fastq_gc_content(fname, ax=axs[1])
            plt.tight_layout()
            w.show_figure(fig)
            i = self.tabs.addTab(w, label )
            self.tabs.setCurrentIndex(i)
        return

    def add_read_lengths(self, progress_callback):
        """get read lengths"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        for i,r in data.iterrows():
            df.loc[i,'reads'] = tools.get_fastq_length(r.filename1)
        return

    def add_mapping_stats(self, progress_callback):
        """get mapping stats for all files and add to table"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]

        def get_stats(x):
            d = tools.samtools_flagstat(x)
            s = pd.Series(d)
            return s

        for i,r in data.iterrows():
            if pd.isnull(r.bam_file):
                print ('no bam file')
                continue
            s = get_stats(r.bam_file)
            df.loc[i,'mapped'] = s['mapped']
            df.loc[i,'total'] = s['total']
        #self.fastq_table.setDataFrame(df)
        return

    def mapping_stats(self, row):
        """Summary of a single fastq file"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        if 'bam_file' not in df.columns:
            print ('No bam file info. Have reads been aligned?')
            return
        if pd.isnull(data.bam_file):
            print ('no bam file')
            return
        d = tools.samtools_flagstat(data.bam_file)
        df = pd.DataFrame(d.items())
        self.info.append(data.bam_file)
        self.info.append(df.to_string())
        self.info.append('-----------------')
        return

    def show_bam_viewer(self, row):
        """Show simple alignment view for a bam file"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = 'aln:'+data['sample']
        if name in self.get_tab_names():
            return
        #text view using samtools tview
        w = widgets.SimpleBamViewer(self)
        w.load_data(data.bam_file, self.ref_genome, self.ref_gb)
        w.redraw(xstart=1)
        i = self.tabs.addTab(w, name )
        self.tabs.setCurrentIndex(i)
        return

    def fastq_quality_report(self):
        """Make fastq quality report as pdf"""

        if self.running == True:
            print ('another process running')
            return
        self.running == True
        df = self.fastq_table.model.df
        out = os.path.join(self.outputdir,'qc_report.pdf')
        if not os.path.exists(out):
            def func(progress_callback):
                tools.pdf_qc_reports(df.filename1, out)
                import webbrowser
                webbrowser.open_new(out)
            self.run_threaded_process(func, self.processing_completed)
        else:
                import webbrowser
                webbrowser.open_new(out)
        return

    def plot_snp_matrix(self):

        mat = pd.read_csv(self.results['snp_dist'],index_col=0)
        bv = widgets.BrowserViewer()
        import toyplot
        min=mat.min().min()
        max=mat.max().max()
        colormap = toyplot.color.brewer.map("BlueGreen", domain_min=min, domain_max=max)
        locator = toyplot.locator.Explicit(range(len(mat)),list(mat.index))
        canvas,axes = toyplot.matrix((mat.values,colormap), llocator=locator, tlocator=locator,
                        label="SNP distance matrix", colorshow=True)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            bv.browser.setHtml(html)

        idx = self.right_tabs.addTab(bv, 'snp dist')
        self.right_tabs.setCurrentIndex(idx)
        return

    def show_browser_tab(self, link, name):

        from PySide2.QtWebEngineWidgets import QWebEngineView
        browser = QWebEngineView()
        browser.setUrl(link)
        idx = self.right_tabs.addTab(browser, name)
        self.right_tabs.setCurrentIndex(idx)
        return

    def check_heterozygosity(self):
        """Plot heterozygosity for each sample"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        samples = list(data['sample'].unique())
        print (samples)
        vcffile = self.results['vcf_file']
        vdf = tools.vcf_to_dataframe(vcffile)
        def het(x):
            if sum(x.AD) == 0:
                return
            return min(x.AD)/sum(x.AD)
        l=int(np.sqrt(len(samples)))
        fig,ax=plt.subplots(l,l+1,figsize=(10,6))
        axs=ax.flat
        i=0
        sites = []
        for s in samples:
            x = vdf[vdf['sample']==s]
            x['het'] = x.apply(het,1)
            x.plot('start','het',kind='scatter',alpha=0.6,ax=axs[i])
            axs[i].set_title(s)
            i+=1
            h = x[x['het']>0.1]
            sites.append(h)

        plt.tight_layout()
        w = widgets.PlotViewer(self)
        w.show_figure(fig)
        i = self.tabs.addTab(w, 'hetero')
        fig.savefig(os.path.join(self.outputdir, 'hetero.png'))
        return

    def snp_typing(self, progress_callback):
        """SNP typing for M.bovis"""

        from . import snp_typing
        df = self.fastq_table.model.df
        #use ALL snp sites including uninformative
        if not 'snp_dist' in self.results:
            print ('You need to create a SNP alignment first')
            return
        snpmat = pd.read_csv(self.results['snp_dist'],sep=' ',index_col=0)
        #print (nucmat)
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        snptable = snp_typing.clade_snps
        res = snp_typing.type_samples(snpmat)
        print (res)
        return

    def spoligotyping(self, progress_callback):
        """Mbovis spo typing tool"""

        self.running == True
        self.opts.applyOptions()
        kwds = self.opts.kwds
        df = self.fastq_table.model.df
        data = self.get_selected()
        if data is None or len(data) == 0:
            return
        res=[]
        cols = ['sample','spotype','sb']
        for i,r in data.iterrows():
            name = r['sample']
            s = tools.get_spoligotype(r.filename1, threads=kwds['threads'])
            sb = tools.get_sb_number(s)
            print (name, s, sb)
            #set new values in place
            df.loc[i,cols] = [name,s,sb]
        return

    def spotyping_completed(self):
        """Typing completed"""

        print("finished")
        self.progressbar.setRange(0,1)
        self.fastq_table.refresh()
        self.running = False
        return

    def rd_analysis(self, progress_callback):
        """Run RD analysis for MTBC species"""

        self.running == True
        self.opts.applyOptions()
        kwds = self.opts.kwds
        from . import rdiff
        rdiff.create_rd_index()

        data = self.get_selected()
        if data is None or len(data) == 0:
            return
        out = os.path.join(self.outputdir,'rd_analysis')
        res = rdiff.run_samples(data, out, threads=kwds['threads'])

        X = rdiff.get_matrix(res, cutoff=0.15)
        X['species'] = X.apply(rdiff.apply_rules,1)
        self.rd_result = X
        return

    def rd_analysis_completed(self):
        """RD analysis completed"""

        print("finished")
        self.progressbar.setRange(0,1)
        if not hasattr(self, 'rd_result'):
            return
        X = self.rd_result
        df = self.fastq_table.model.df

        right = X[['species']]
        #df = df.merge(X[['species']],left_on='sample',right_index=True,how='left')
        df = df.set_index('sample').join(right,rsuffix='_x')
        if 'species_x' in df.columns:
            df['species'] = df.species.fillna(df.species_x)
            df = df.drop('species_x', axis=1)
        df.reset_index(inplace=True)
        self.fastq_table.model.df = df
        self.fastq_table.refresh()
        #add plot
        fig,ax = plt.subplots(1,1)
        plotting.heatmap(X.set_index('species',append=True), cmap='cubehelix',ax=ax)
        w = widgets.PlotViewer(self)
        w.show_figure(fig)
        i = self.tabs.addTab(w, 'RD')
        self.running = False
        return

    def get_selected(self):
        """Get selected rows of fastq table"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        if len(rows) == 0:
            print ('no samples selected')
            return
        data = df.iloc[rows]
        return data

    def zoom_in(self):

        w = self.fastq_table
        w.zoomIn()
        self.info.zoomIn()
        return

    def zoom_out(self):

        w = self.fastq_table
        w.zoomOut()
        self.info.zoomOut()
        return

    def show_info(self, msg, color=None):

        if color != None:
            #self.info.appendHtml("<p style=\"color:red\">" + msg + "</p>")
            self.info.append("<font color=%s>%s</font>" %(color,msg))
        else:
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
        self.save_settings()
        event.accept()
        return

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def online_documentation(self,event=None):
        """Open the online documentation"""

        #import webbrowser
        link='https://github.com/dmnfarrell/snipgenie'
        #webbrowser.open(link,autoraise=1)
        from PySide2.QtWebEngineWidgets import QWebEngineView
        browser = QWebEngineView()
        browser.setUrl(link)
        idx = self.right_tabs.addTab(browser, 'help')
        self.right_tabs.setCurrentIndex(idx)
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

        text='snipgenie GUI\n'\
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

class StdoutRedirect(QObject):
    printOccur = Signal(str, str, name="print")

    def __init__(self, *param):
        QObject.__init__(self, None)
        self.daemon = True
        self.sysstdout = sys.stdout.write
        self.sysstderr = sys.stderr.write

    def stop(self):
        sys.stdout.write = self.sysstdout
        sys.stderr.write = self.sysstderr

    def start(self):
        sys.stdout.write = self.write
        sys.stderr.write = lambda msg : self.write(msg, color="red")

    def write(self, s, color="black"):
        sys.stdout.flush()
        self.printOccur.emit(s, color)

class AppOptions(widgets.BaseOptions):
    """Class to provide a dialog for global plot options"""

    def __init__(self, parent=None):
        """Setup variables"""
        self.parent = parent
        self.kwds = {}
        genomes = []
        aligners = ['bwa','subread','minimap2']
        platforms = ['illumina','ont']
        separators = ['_','-','|',';','~']
        cpus = os.cpu_count()
        self.groups = {'general':['threads','labelsep','overwrite'],
                        'trimming':['quality'],
                        'aligners':['aligner','platform'],
                        'variant calling':['filters','proximity'],
                        #'blast':['db','identity','coverage']
                       }
        self.opts = {'threads':{'type':'spinbox','default':4,'range':(1,cpus)},
                    'overwrite':{'type':'checkbox','default':False},
                    'labelsep':{'type':'combobox','default':'_',
                    'items':separators,'label':'label sep','editable':True},
                    'aligner':{'type':'combobox','default':'bwa',
                    'items':aligners,'label':'aligner'},
                    'platform':{'type':'combobox','default':'illumina',
                    'items':platforms,'label':'platform'},
                    'filters':{'type':'entry','default':app.default_filter},
                    'proximity': {'type':'checkbox','default':0},
                    'quality':{'type':'spinbox','default':30}
                    #'db':{'type':'combobox','default':'card',
                    #'items':[],'label':'database'},
                    #'identity':{'type':'entry','default':90},
                    #'coverage':{'type':'entry','default':50},
                    }
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snipgenie gui tool')
    parser.add_argument("-f", "--fasta", dest="filenames",default=[],
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-p", "--proj", dest="project",default=None,
                        help="load .snipgenie project file", metavar="FILE")
    args = vars(parser.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = App(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
