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
import matplotlib as mpl
from . import core, tools, aligners, app, widgets, tables, plotting, trees

homepath = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
logoimg = os.path.join(module_path, 'logo.svg')
stylepath = os.path.join(module_path, 'styles')
iconpath = os.path.join(module_path, 'icons')
pluginiconpath = os.path.join(module_path, 'plugins', 'icons')
settingspath = os.path.join(homepath, '.config','snipgenie')

dockstyle = '''
    QDockWidget {
        max-width:240px;
    }
    QDockWidget::title {
        background-color: #ccd8f0;
    }
    QScrollBar:vertical {
         width: 15px;
         margin: 1px 0 1px 0;
     }
    QScrollBar::handle:vertical {
         min-height: 20px;
     }
'''

widgetstyle = '''
    QWidget {
        font-size: 12px;
        max-width: 220px;
    }
    QLabel {
        min-width: 60px;
        width:80px;
    }
    QPlainTextEdit {
        max-height: 100px;
        min-width: 100px;
    }
    QScrollBar:vertical {
         width: 15px;
     }
    QComboBox {
        combobox-popup: 0;
        max-height: 30px;
        max-width: 100px;
    }
    QListView::item:selected {
        min-width: 300px;}
'''

#fix to QtWebEngine display on linux
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--no-sandbox"

class Communicate(QObject):
    newproj = Signal()

class App(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, filenames=[], project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.comms = Communicate()
        self.setWindowTitle("snipgenie")

        self.setWindowIcon(QIcon(logoimg))
        self.create_menu()
        self.main = QSplitter()

        screen_resolution = QGuiApplication.primaryScreen().availableGeometry()
        width, height = screen_resolution.width()*0.7, screen_resolution.height()*.7
        if screen_resolution.width()>1280:
            self.setGeometry(QtCore.QRect(200, 200, int(width), int(height)))
        else:
            self.showMaximized()
        self.setMinimumSize(400,300)

        self.running = False
        self.recent_files = ['']
        self.scratch_items = {}
        self.opentables = {}
        self.openplugins = {}
        self.plugindata = {}

        self.main.setFocus()
        self.setCentralWidget(self.main)

        self.create_tool_bar()
        self.setup_gui()
        self.load_settings()
        self.show_recent_files()
        self.start_logging()

        self.new_project()

        if platform.system() == 'Windows':
            app.fetch_binaries()
        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
        mpl.style.use('bmh')
        self.discover_plugins()
        self.redirect_stdout()
        return

    def start_logging(self):
        """Error logging"""

        import logging
        if platform.system() == 'Windows':
            path = QtCore.QStandardPaths.writableLocation(QtCore.QStandardPaths.ConfigLocation)
            if not os.path.exists(path):
                os.makedirs(path)
        else:
            path = os.path.dirname(self.settings.fileName())
        self.logfile = os.path.join(path, 'error.log')
        logging.basicConfig(filename=self.logfile,format='%(asctime)s %(message)s')
        return

    def redirect_stdout(self):
        """redirect stdout"""
        self._stdout = StdoutRedirect()
        self._stdout.start()
        self._stdout.printOccur.connect(lambda x : self.info.insert(x))
        return

    def load_settings(self):
        """Load GUI settings"""

        s = self.settings = QtCore.QSettings('snipgenie','default')
        try:
            winsize = s.value('window_size')
            if winsize != None:
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

    def apply_settings(self):
        """Apply settings to GUI when changed"""

        self.setIconSize(QtCore.QSize(core.ICONSIZE, core.ICONSIZE))
        for i in self.opentables:
            table = self.opentables[i]
            #table.toolbar.setIconSize(QtCore.QSize(core.ICONSIZE, core.ICONSIZE))
        import matplotlib as mpl
        mpl.rcParams['savefig.dpi'] = core.DPI
        #self.setTheme(self.theme)
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
                    'action': lambda: self.run_threaded_process(self.variant_calling, self.calling_completed),
                    'file':'call-variants'},
                 'SNP Dist Matrix': {
                    'action':  self.show_snpdist,
                    'file':'snp-dist'},
                 'SNP Viewer': {
                    'action':  self.snp_viewer,
                    'file':'snp-viewer'},
                 'Show Phylogeny': {
                    'action': self.show_phylogeny,
                    'file':'phylogeny'},
                 'Scratchpad': {'action':self.show_scratchpad,'file':'scratchpad'},
                 'Quit': {'action':self.quit,'file':'application-exit'}
                }

        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
        for i in items:
            if 'file' in items[i]:
                iconfile = os.path.join(iconpath,items[i]['file'])
                icon = QIcon(iconfile)
            else:
                icon = QIcon.fromTheme(items[i]['icon'])
            btn = QAction(icon, i, self)
            btn.triggered.connect(items[i]['action'])
            #btn.setCheckable(True)
            toolbar.addAction(btn)
        return

    def add_dock(self, widget, name):
        """Add a dock widget"""

        dock = QDockWidget(name)
        dock.setStyleSheet(dockstyle)
        area = QScrollArea()
        area.setWidgetResizable(True)
        dock.setWidget(area)
        area.setWidget(widget)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        self.docks['options'] = dock
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.docks = {}
        #self.m = QSplitter(self.main)
        self.m = self.main

        style = '''
        QWidget {
            font-size: 12px;
            max-height: 160px;
            }
        '''
        dialog = QWidget()
        dialog.setStyleSheet(style)
        #dialog.setFixedSize(200, 200)
        l = QVBoxLayout(dialog)
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

        self.add_dock(dialog, 'genome')
        #create option widgets
        self.opts = AppOptions(parent=self.m)
        dialog = self.opts.showDialog(self, wrap=1, section_wrap=1, style=widgetstyle)
        self.add_dock(dialog, 'options')

        #add dock menu items
        for name in ['options']:
            action = self.docks[name].toggleViewAction()
            self.dock_menu.addAction(action)
            action.setCheckable(True)

        #general plot window
        #self.plotview = widgets.PlotViewer(self)

        center = QSplitter(Qt.Vertical)
        self.m.addWidget(center)

        self.fastq_table = tables.SampleTable(self, dataframe=pd.DataFrame, app=self)
        self.table_widget = tables.DataFrameWidget(parent=center, table=self.fastq_table,
                            toolbar=True)
        #self.fastq_table = self.table_widget.table
        center.addWidget(self.table_widget)
        self.opentables['main'] = self.fastq_table

        self.tabs = QTabWidget(center)
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        center.addWidget(self.tabs)
        center.setSizes((100,100))

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

        #self.right_tabs.addTab(self.plotview, 'plots')

        self.m.setSizes([200,150])
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
        return

    @Slot(int)
    def close_tab(self, index):
        """Close current tab"""

        #index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        self.tabs.removeTab(index)
        return

    @Slot(int)
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

        self.edit_menu = QMenu('Edit', self)
        self.menuBar().addMenu(self.edit_menu)
        #self.edit_menu.addAction(icon, 'Find/Replace', self.findReplace,
        #        QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        icon = QIcon(os.path.join(iconpath,'preferences-system.png'))
        self.edit_menu.addAction(icon, 'Preferences', self.preferences)

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

        self.analysis_menu = QMenu('Workflow', self)
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
        #self.analysis_menu.addAction('Create SNP alignment',
        #    lambda: self.run_threaded_process(self.snp_alignment, self.snp_align_completed))
        icon = QIcon(os.path.join(iconpath,'phylogeny.png'))
        self.analysis_menu.addAction(icon, 'Build Phylogeny',
            lambda: self.run_threaded_process(self.make_phylo_tree, self.phylogeny_completed))
        #self.analysis_menu.addSeparator()
        #self.analysis_menu.addAction('Run Workflow', self.run)

        self.tools_menu = QMenu('Tools', self)
        self.menuBar().addMenu(self.tools_menu)

        self.tools_menu.addAction('Get Read Stats',
            lambda: self.run_threaded_process(self.add_read_lengths, self.processing_completed))
        self.tools_menu.addAction('Get Mapping Stats',
            lambda: self.run_threaded_process(self.add_mapping_stats, self.processing_completed))
        self.tools_menu.addAction('Get Depth/Coverage',
            lambda: self.run_threaded_process(self.add_mean_depth, self.processing_completed))
        self.tools_menu.addAction('Missing Sites',
            lambda: self.run_threaded_process(self.missing_sites, self.processing_completed))
        self.tools_menu.addAction('Mean GC content',
            lambda: self.run_threaded_process(self.add_gc_mean, self.processing_completed))
        self.tools_menu.addAction('Fastq Qualities Report', self.fastq_quality_report)
        self.tools_menu.addAction('Consensus Sequences', self.get_consensus_sequences)

        self.tools_menu.addAction('Show Annotation', self.show_ref_annotation)
        self.tools_menu.addAction('SNP Dist Matrix', self.show_snpdist)
        self.tools_menu.addAction('SNP Viewer', self.snp_viewer)
        self.tools_menu.addAction('CSQ Viewer', self.csq_viewer)
        self.tools_menu.addAction('VCF Viewer', self.vcf_viewer)
        #self.tools_menu.addAction('Map View', self.show_map)
        self.tools_menu.addAction('Tree Viewer', self.tree_viewer)
        #self.tools_menu.addSeparator()
        #self.tools_menu.addAction('Check Heterozygosity', self.check_heterozygosity)
        #self.tools_menu.addAction('RD Analysis (MTBC)',
        #    lambda: self.run_threaded_process(self.rd_analysis, self.rd_analysis_completed))
        #self.tools_menu.addAction('M.bovis SNP typing',
        #    lambda: self.run_threaded_process(self.snp_typing, self.processing_completed))

        self.settings_menu = QMenu('Settings', self)
        self.menuBar().addMenu(self.settings_menu)
        self.settings_menu.addAction('Set Output Folder', self.set_output_folder)
        self.settings_menu.addAction('Set Reference Sequence', self.set_reference)
        self.settings_menu.addAction('Set Annnotation (genbank)', self.set_annotation)
        #self.settings_menu.addAction('Set Filters', self.set_filters)
        self.settings_menu.addAction('Add Mask File', self.add_mask)
        self.settings_menu.addAction('Add Sample Meta Data', self.merge_meta_data)
        self.settings_menu.addAction('Clean Up Files', self.clean_up)

        self.presets_menu = QMenu('Preset Genomes', self)
        self.menuBar().addMenu(self.presets_menu)
        self.load_presets_menu()

        #self.tools_menu.setEnabled(False)

        self.dock_menu = QMenu('Docks', self)
        self.menuBar().addMenu(self.dock_menu)

        self.plugin_menu = QMenu('Plugins', self)
        self.menuBar().addMenu(self.plugin_menu)

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('View Error Log', self.show_error_log)
        self.help_menu.addAction('Help', self.online_documentation)
        self.help_menu.addAction('Online BLAST', self.show_blast_url)
        self.help_menu.addAction('Nucleotide DB search', self.show_nucldb_url)
        self.help_menu.addAction('About', self.about)

    def load_presets_menu(self,ask=True):
        """Add preset genomes to menu"""

        genomes = app.preset_genomes
        for name in genomes:
            seqname = genomes[name]['sequence']
            gbfile = genomes[name]['gb']
            if 'mask' in genomes[name]:
                mask = genomes[name]['mask']
            else:
                mask = None
            '''def func(seqname,gbfile,mask,ask):

                if ask == True:
                    msg = 'This will replace the current reference with a preset.'
                    reply = QMessageBox.question(self, 'Warning!', msg,
                                                    QMessageBox.No | QMessageBox.Yes )
                    if reply == QMessageBox.No:
                        return
                self.set_reference(seqname, ask=False)
                self.set_annotation(gbfile)
                self.set_mask(mask)'''
            func = self.load_preset_genome
            receiver = lambda seqname=seqname, gb=gbfile, mask=mask, ask=ask: func(seqname, gb, mask, ask)
            self.presets_menu.addAction('%s' %name, receiver)
        return

    def refresh(self):
        """Refresh all tables"""

        for i in self.opentables:
            w = self.opentables[i]
            w.font = core.FONT
            w.fontsize = core.FONTSIZE
            w.refresh()
        return

    def load_preset_genome(self, seqname, gbfile, mask, ask):

        if ask == True:
            msg = 'This will replace the current reference with a preset.'
            reply = QMessageBox.question(self, 'Warning!', msg,
                                            QMessageBox.No | QMessageBox.Yes )
            if reply == QMessageBox.No:
                return
        self.set_reference(seqname, ask=False)
        self.set_annotation(gbfile)
        self.set_mask(mask)
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
        data['opentables'] = list(self.opentables.keys())
        self.opts.applyOptions()
        data['options'] = self.opts.kwds
        data['scratch_items'] = self.scratch_items
        self.projectlabel.setText(filename)
        data['plugins'] = self.save_plugin_data()
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

        if self.running == True:
            QMessageBox.information(self, "Process is Running",
                            'Wait until process is finished.')
            return False
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
        self.scratch_items = {}
        self.proj_file = None
        self.fastq_table.setDataFrame(pd.DataFrame({'name':[]}))
        self.fastq_table.refresh()
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
        self.clear_plugins()
        self.comms.newproj.emit()
        return True

    def setup_paths(self):
        """Set paths to important files in proj folder"""

        self.core_fasta = os.path.join(self.outputdir, 'core.fa')
        self.dist_matrix = os.path.join(self.outputdir, 'snpdist.csv')
        self.snp_matrix = os.path.join(self.outputdir, 'core.txt')
        self.csq_matrix = os.path.join(self.outputdir, 'csq.matrix')
        self.treefile = os.path.join(self.outputdir,'tree.newick')
        return

    def load_project(self, filename=None):
        """Load project"""

        closed = self.new_project()
        if closed == False:
            return
        data = pickle.load(open(filename,'rb'))
        keys = ['sheets','outputdir','results','ref_genome','ref_gb','mask_file']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        if 'options' in data:
            self.opts.updateWidgets(data['options'])
        if 'plugins' in data:
            self.plugindata = data['plugins']
        else:
            self.plugindata = {}
        ft = self.fastq_table
        ft.setDataFrame(data['inputs'])
        ft.resizeColumns()
        ft.refresh()

        self.check_files()
        self.update_ref_genome()
        self.update_mask()
        self.proj_file = filename
        self.update_labels()
        self.setup_paths()
        #self.show_snpdist()
        if 'scratch_items' in data:
            self.scratch_items = data['scratch_items']
        #if 'opentables' in data:
            #for key in data['opentables']:
                #print (key)
        #load tree view
        if 'treeviewer' in data.keys():
            self.tree_viewer()
            self.treeviewer.loadData(data['treeviewer'])
        self.right_tabs.setCurrentIndex(0)
        self.add_recent_file(filename)
        self.check_fastq_table()
        return

    def update_labels(self):
        self.projectlabel.setText(self.proj_file)
        self.outdirLabel.setText(self.outputdir)

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

    def check_files(self):
        """Check input files exist"""

        ft = self.fastq_table
        df = ft.model.df
        for i,r in df.iterrows():
            if not os.path.exists(r.filename1):
                print ('%s missing file' %i)
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
        cols = ['bam_file','mapped','reads']
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
        app.write_samples(df[['sample']], self.outputdir)
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
        print (kwds)
        df = self.fastq_table.model.df
        new = app.get_samples(filenames, sep=kwds['labelsep'])
        #pivoted
        new = app.get_pivoted_samples(new)
        print ('getting fastq read lengths..')
        new['read_length'] = new.filename1.apply(lambda x: tools.get_fastq_info(x))

        if len(df)>0:
            new = pd.concat([df,new],sort=False).reset_index(drop=True)
            new = new.drop_duplicates('filename1')

        self.fastq_table.setDataFrame(new)
        self.fastq_table.resizeColumns()
        #print (new)
        app.write_samples(new[['sample']], self.outputdir)
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

    def set_reference(self, filename=None, ask=True):
        """Reset the reference sequence"""

        msg = "This will change the reference genome. You will need to re-run any previous alignments. Are you sure?"
        if ask == True:
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
        if not filename:
            return
        self.ref_gb = filename
        #put annotation in a dataframe
        self.annot = tools.genbank_to_dataframe(self.ref_gb)
        self.update_ref_genome()
        return

    def set_mask(self, filename):

        self.mask_file = filename
        self.update_mask()

    def add_mask(self, filename=None):
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
                #ensure sample col object type when we import
                df = pd.read_csv(results_file, dtype={'sample':'object'})
                self.fastq_table.model.df = df
                self.fastq_table.refresh()
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
        #print (rows)
        new = df.iloc[rows]
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
        self.snp_alignment()

        return

    def calling_completed(self):

        self.processing_completed()
        self.show_snpdist()
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

        if not os.path.exists(self.dist_matrix):
            return
        mat = pd.read_csv(self.dist_matrix,index_col=0)
        if 'SNP dist' in self.get_tabs():
            self.tabs.removeTab(0)
        table = tables.DistMatrixTable(self, app=self, dataframe=mat)
        w = tables.DataFrameWidget(self.tabs, table=table,
                            toolbar=True)
        i = self.tabs.addTab(w, 'SNP dist')
        self.tabs.setCurrentIndex(i)
        self.opentables['SNP dist'] = table
        return

    def snp_alignment(self, progress_callback=None):
        """Make snp matrix from variant positions"""

        self.opts.applyOptions()
        kwds = self.opts.kwds
        vcf_file = self.results['vcf_file']
        print('Making SNP alignment')
        result, smat = tools.core_alignment_from_vcf(vcf_file)
        #print (result)
        outfasta = os.path.join(self.outputdir, 'core.fa')
        SeqIO.write(result, outfasta, 'fasta')
        smat.to_csv(os.path.join(self.outputdir,'core.txt'), sep=' ')
        from Bio import AlignIO
        aln = AlignIO.read(outfasta, 'fasta')
        snp_dist = tools.snp_dist_matrix(aln)
        snp_dist.to_csv(self.dist_matrix, sep=',')
        return

    def missing_sites(self, progress_callback=None):
        """Find missing sites in each sample - useful for quality control"""

        vcf_file = os.path.join(self.outputdir, 'snps.vcf.gz')
        snprecs, smat = tools.core_alignment_from_vcf(vcf_file, missing=True)
        #outfasta = os.path.join(self.outdir, 'core.fa')
        #SeqIO.write(snprecs, outfasta, 'fasta')
        #write out sites matrix as txt file
        smat.to_csv(os.path.join(self.outputdir,'core_missing.txt'), sep=' ')
        missing = smat[smat=='N'].count().sort_values()

        #merge with samples table - bad!
        df = self.fastq_table.model.df
        df = df.merge(missing.rename('missing_sites'),how='left',
                        right_index=True,left_on='sample')
        self.fastq_table.model.df = df
        self.fastq_table.refresh()
        return

    def get_tab_indices(self, tab_widget, tab_name):
        return [index for index in range(tab_widget.count())
            if tab_name == tab_widget.tabText(index)]

    def snp_viewer(self):
        """Show SNP table - output of core.txt"""

        file = os.path.join(self.outputdir, 'core.txt')
        mat = pd.read_csv(file, sep=' ',index_col=0).sort_index()
        mat = mat.T
        if 'SNP' in self.get_tabs():
            index = self.get_tab_indices(self.tabs, 'SNP')
            print (index)
            self.tabs.removeTab(index[0])
        table = tables.SNPTable(self.tabs, app=self, dataframe=mat)
        idx = self.tabs.addTab(table, 'SNP')
        self.tabs.setCurrentIndex(idx)
        self.opentables['SNP'] = table
        return

    def csq_viewer(self):
        """Show CSQ table - output of bcftools csq"""

        if not os.path.exists(self.csq_matrix):
            return
        mat = pd.read_csv(self.csq_matrix)
        if 'CSQ' in self.get_tabs():
            self.tabs.removeTab(0)
        table = tables.CSQTable(self.tabs, app=self, dataframe=mat)
        idx = self.tabs.addTab(table, 'CSQ')
        self.tabs.setCurrentIndex(idx)
        self.opentables['CSQ'] = table
        return

    def vcf_viewer(self):
        """Show VCF table"""

        if 'VCF' in self.get_tabs():
            self.tabs.removeTab(0)
        vcf_file = os.path.join(self.outputdir, 'snps.vcf.gz')
        df = tools.vcf_to_dataframe(vcf_file).set_index('sample')
        table = tables.VCFTable(self.tabs, app=self, dataframe=df)
        idx = self.tabs.addTab(table, 'VCF')
        self.tabs.setCurrentIndex(idx)
        self.opentables['VCF'] = table
        return

    def make_phylo_tree(self, progress_callback=None, method='raxml'):
        """Make phylogenetic tree"""

        corefasta = os.path.join(self.outputdir, 'core.fa')
        bootstraps = 100
        if method == 'raxml':
            treefile = trees.run_RAXML(corefasta, bootstraps=bootstraps, outpath=self.outputdir)
        elif method == 'fasttree':
            treefile = trees.run_fasttree(corefasta, bootstraps=bootstraps, outpath=self.outputdir)
        outfile = os.path.join(self.outputdir,'tree.newick')
        snpmat = pd.read_csv(self.snp_matrix, sep=' ',index_col=0)
        ls = len(snpmat)
        trees.convert_branch_lengths(treefile, outfile, ls)
        self.treefile = outfile
        return

    def phylogeny_completed(self):

        self.processing_completed()
        self.show_phylogeny()

    def show_phylogeny(self):
        """Show current tree"""

        self.tree_viewer()
        filename = os.path.join(self.outputdir,'tree.newick')
        self.treeviewer.load_tree(filename)
        self.treeviewer.update()
        return

    def tree_viewer(self):
        """Show tree viewer"""

        from . import phylo
        if not hasattr(self, 'treeviewer'):
            self.treeviewer = phylo.TreeViewer(self)
        if not 'phylogeny' in self.get_right_tabs():
            idx = self.right_tabs.addTab(self.treeviewer, 'phylogeny')
            self.right_tabs.setCurrentIndex(idx)
        return

    def show_map(self):

        from . import gis
        if not hasattr(self, 'gisviewer'):
            self.gisviewer = gis.GISViewer()
        if not 'map' in self.get_right_tabs():
            idx = self.right_tabs.addTab(self.gisviewer, 'map')
            self.right_tabs.setCurrentIndex(idx)
        return

    def get_tabs(self):

        n=[]
        for i in range(self.tabs.count()):
            n.append(self.tabs.tabText(i))
        return n

    def get_right_tabs(self):

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
        self.running = True
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
        #print (data)
        #print ()
        w = widgets.TableViewer(self, pd.DataFrame(data))
        i = self.right_tabs.addTab(w, data['sample'])
        self.right_tabs.setCurrentIndex(i)
        return

    def quality_summary(self, row):
        """Summary of a single sample, both files"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        rl = tools.get_fastq_read_lengths(data.filename1)
        if 'filename2' in df.columns:
            colnames = zip([data.name1, data.name2], [data.filename1, data.filename2])
        else:
            colnames = [[data.name1, data.filename1]]
        for name,fname in colnames:
            label = 'quality:'+name
            if label in self.get_tab_names():
                return
            w = widgets.PlotViewer(self)
            fig,ax = plt.subplots(3,1, figsize=(7,5), dpi=100,
                        facecolor=(1,1,1), edgecolor=(0,0,0))
            axs=ax.flat
            if not os.path.exists(fname):
                self.show_info('This file is missing.')
                return
            if rl.mean()<800:
                tools.plot_fastq_qualities(fname, ax=axs[0])
            tools.plot_fastq_gc_content(fname, ax=axs[1])

            c = app.blast_contaminants(fname)
            #c.perc_hits.plot(kind='pie',ax=axs[2])
            ax=axs[2]
            c.perc_hits.plot(kind='bar',ax=ax)
            labels = ax.get_xticklabels()
            labels = [label.get_text().replace('_', '\n') for label in labels ]
            ax.set_xticklabels(labels, rotation=0)
            ax.set_ylabel('% total')
            ax.set_title('contaminant check')

            fig.suptitle('Qualities: %s' %name, fontsize=18)
            plt.tight_layout()
            w.set_figure(fig)
            i = self.right_tabs.addTab(w, label)
            self.right_tabs.setCurrentIndex(i)
        return

    def read_distributon(self, row):
        """get read length distribution"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        name=data['sample']
        rl = tools.get_fastq_read_lengths(data.filename1)

        w = widgets.PlotViewer(self)
        fig,ax = plt.subplots(1,1, figsize=(7,5), dpi=65)
        rl.hist(bins=20,ax=ax)
        ax.set_title(name)
        w.show_figure(fig)
        label = 'readlengths:'+name
        i = self.tabs.addTab(w, label )
        self.tabs.setCurrentIndex(i)
        return

    def add_read_lengths(self, progress_callback):
        """Get read lengths"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        for i,r in data.iterrows():
            total = tools.get_fastq_size(r.filename1)
            df.loc[i,'reads'] = total
        return

    def add_gc_mean(self, progress_callback):
        """Get mean GC to indicate contamination"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        for i,r in data.iterrows():
            df.loc[i,'meanGC'] = tools.get_gc(r.filename1, limit=5e4).mean().round(2)
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
            df.loc[i,'mapped'] = s['primary']
            total = df.loc[i,'reads']
            if pd.isna(total):
                total = tools.get_fastq_size(r.filename1)
                df.loc[i,'reads'] = total
            df.loc[i,'perc_mapped'] = round(s['primary']/(total*2)*100,2)
            #print (s['mapped'],total)
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

    def add_mean_depth(self, progress_callback):
        """find mean depth for bam file"""

        df = self.fastq_table.model.df
        if 'bam_file' not in df.columns:
            print ('No bam file info. Have reads been aligned?')
            return
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        cols = ['coverage','meandepth','meanbaseq','meanmapq']
        for i,r in data.iterrows():
            if pd.isnull(r.bam_file):
                continue
            #df.loc[i,'depth'] = tools.get_mean_depth(r.bam_file)
            c = tools.samtools_coverage(r.bam_file)
            df.loc[i,cols] = c[cols]
        return

    def get_consensus_sequences(self):
        """Consensus sequences for all samples"""

        vcf_file = os.path.join(self.outputdir, 'snps.vcf.gz')
        outfile = os.path.join(self.outputdir, 'consensus.fa')
        tools.bcftools_consensus(vcf_file, self.ref_genome, outfile)
        return

    def show_bam_viewer(self, row):
        """Show simple alignment view for a bam file"""

        if self.ref_genome == None:
            self.info.append('no reference genome set!')
            return
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

    def check_contamination(self):
        """Blast to common contaminant sequences"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        name = data['sample']
        c = app.blast_contaminants(data.filename1, limit=10000, pident=95)
        if len(c) == 0:
            print ('no contaminants in DB found')
            return
        w = widgets.PlotViewer(self)
        #fig,ax = plt.subplots(1,1, figsize=(7,5), dpi=120)
        ax = w.ax
        c.perc_hits.plot(kind='barh',ax=ax)
        ax.set_xlabel('% total')
        ax.set_title('sequence contaminants - %s' %name)
        plt.tight_layout()
        #w.set_figure(fig)

        i = self.tabs.addTab(w, 'contam:'+name)
        self.tabs.setCurrentIndex(i)
        return

    def get_fasta_reads(self):
        """Get a sample of reads for blasting"""

        df = self.fastq_table.model.df
        row = self.fastq_table.getSelectedRows()[0]
        data = df.iloc[row]
        seqs = tools.fastq_to_rec(data.filename1, size=50)
        fastafmt = '\n'.join([seq.format("fasta") for seq in seqs])
        name = 'seqs:'+data['sample']
        w = widgets.Editor(self, readOnly=True, fontsize=11)
        w.append(fastafmt)
        i = self.right_tabs.addTab(w, name)
        self.right_tabs.setCurrentIndex(i)
        return

    def show_blast_url(self):

        link = 'https://blast.ncbi.nlm.nih.gov/'
        self.show_browser_tab(link, 'NCBI Blast')
        return

    def show_nucldb_url(self):

        link = 'https://www.ncbi.nlm.nih.gov/nuccore/'
        self.show_browser_tab(link, 'NCBI Blast')
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

    def plot_dist_matrix(self):

        mat = pd.read_csv(self.dist_matrix,index_col=0)
        #fig,ax = plt.subplots(1,1,figsize=(8,8))
        #plotting.heatmap(mat, cmap='coolwarm', ax=ax)
        import seaborn as sns
        cm = sns.clustermap(mat, cmap='coolwarm', annot=True, fmt='.0f')
        fig = cm.fig
        plt.tight_layout()
        w = widgets.PlotViewer(self)
        w.set_figure(fig)
        idx = self.right_tabs.addTab(w, 'SNP dist')
        self.right_tabs.setCurrentIndex(idx)
        return

    def show_scratchpad(self):

        if not hasattr(self, 'scratchpad'):
            self.scratchpad = widgets.ScratchPad()
            try:
                self.scratchpad.resize(self.settings.value('scratchpad_size'))
            except:
                pass
        self.scratchpad.update(self.scratch_items)
        self.scratchpad.show()
        self.scratchpad.activateWindow()
        return

    def show_browser_tab(self, link, name):
        """Show web page in a tab"""

        #from PySide2.QtWebEngineWidgets import QWebEngineView
        #browser = QWebEngineView()
        bw = widgets.BrowserViewer(self)
        url = QUrl.fromUserInput(link)
        bw.browser.setUrl(url)

        idx = self.right_tabs.addTab(bw, name)
        self.right_tabs.setCurrentIndex(idx)
        return

    '''def check_heterozygosity(self):
        """Plot heterozygosity for each sample"""

        df = self.fastq_table.model.df
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        samples = list(data['sample'].unique())
        #print (samples)
        vcffile = self.results['vcf_file']
        vdf = tools.vcf_to_dataframe(vcffile)
        def het(x):
            if sum(x.AD) == 0:
                return
            return min(x.AD)/sum(x.AD)
        l = int(np.sqrt(len(samples)))

        fig,ax = plt.subplots(l,l,figsize=(10,6))
        if l==1:
            axs=[ax]
        else:
            axs=ax.flat
        i=0
        sites = []
        for s in samples:
            x = vdf[vdf['sample']==s].copy()
            x['het'] = x.apply(het,1)
            x.plot('start','het',kind='scatter',alpha=0.6,ax=axs[i])
            axs[i].set_title(s)
            i+=1
            h = x[x['het']>0.1]
            sites.append(h)

        plt.tight_layout()
        w = widgets.PlotViewer(self)
        w.set_figure(fig)
        i = self.tabs.addTab(w, 'hetero')
        #fig.savefig(os.path.join(self.outputdir, 'hetero.png'))
        return'''

    '''def snp_typing(self, progress_callback):
        """SNP typing for M.bovis"""

        from . import snp_typing
        df = self.fastq_table.model.df
        #use ALL snp sites including uninformative
        if not os.path.exists(self.snp_matrix):
            print ('You need to create a SNP alignment first')
            return
        snpmat = pd.read_csv(self.snp_matrix, sep=' ',index_col=0)
        #print (nucmat)
        rows = self.fastq_table.getSelectedRows()
        data = df.iloc[rows]
        snptable = snp_typing.clade_snps
        res = snp_typing.type_samples(snpmat)
        print (res)
        return'''

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
        out = os.path.join(self.outputdir,'rd_aligned')
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

        for i in self.opentables:
            w=self.opentables[i]
            w.zoomIn()
        self.info.zoomIn()
        return

    def zoom_out(self):

        for i in self.opentables:
            w=self.opentables[i]
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
        """Close main window"""

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
        QApplication.closeAllWindows()
        return

    def _check_snap(self):
        if os.environ.has_key('SNAP_USER_COMMON'):
            print ('running inside snap')
            return True
        return False

    def show_error_log(self):
        """Show log file contents"""

        f = open(self.logfile,'r')
        s = ''.join(f.readlines())
        dlg = widgets.TextViewer(self, s, title='Log', width=800, height=400)
        dlg.exec_()
        return

    def discover_plugins(self):
        """Discover available plugins"""

        from . import plugin
        default = os.path.join(module_path, 'plugins')
        other = os.path.join(settingspath, 'plugins')
        paths = [default,other]
        failed = plugin.init_plugin_system(paths)
        self.update_plugin_menu()
        return

    def load_plugin(self, plugin):
        """Instantiate the plugin and show widget or run it"""

        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        if not hasattr(self, 'openplugins'):
            self.openplugins = {}
        openplugins = self.openplugins
        c = plugin.capabilities

        if plugin.name in openplugins:
            p = openplugins[plugin.name]
            self.docks[plugin.name].show()
        elif 'gui' in plugin.capabilities:
            try:
                p = plugin(parent=self)
                #track which plugin is running
                openplugins[plugin.name] = p
                #load data if any saved in project
                #print (self.plugindata)
                if plugin.name in self.plugindata:
                    p.load_data(self.plugindata[plugin.name])
            except Exception as e:
                QMessageBox.information(self, "Plugin error", str(e))
                print(traceback.format_exc())
                return
            #plugin should be added as a dock widget
            self.show_plugin(p)
        else:
            #no widgets, just run the plugins run method
            p = plugin(parent=self)
            p.run()
        return

    def show_plugin(self, plugin):
        """Show plugin in dock or as window"""

        c = plugin.capabilities
        if 'docked' in c:
            self.add_plugin_dock(plugin)
        else:
            plugin.main.show()
            plugin.main.setWindowTitle(plugin.name)
        self.comms.newproj.connect(plugin.project_closed)
        return

    def add_plugin_dock(self, plugin):
        """Add plugin as dock widget"""

        dockstyle = '''
            QDockWidget::title {
                background-color: #EDC6B3;
            }
        '''
        dock = QDockWidget(plugin.name)
        dock.setStyleSheet(dockstyle)
        area = QScrollArea()
        area.setWidgetResizable(True)
        dock.setWidget(area)
        area.setWidget(plugin.main)
        if plugin.side == 'left':
            self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        else:
            self.addDockWidget(QtCore.Qt.RightDockWidgetArea, dock)
        if plugin.side == 'floating':
            dock.setFloating(True)
            #dock.setGeometry(100, 0, 200, 30)
        self.docks[plugin.name] = dock
        return

    def update_plugin_menu(self):
        """Update plugins"""

        from . import plugin
        plgmenu = self.plugin_menu
        #for plg in plugin.get_plugins_classes('gui'):
        for plg in plugin.Plugin.__subclasses__():
            #print (plg)
            def func(p, **kwargs):
                def new():
                   self.load_plugin(p)
                return new
            if hasattr(plg,'iconfile'):
                icon = QIcon(os.path.join(pluginiconpath,plg.iconfile))
                plgmenu.addAction(icon, plg.menuentry, func(plg))
            else:
                plgmenu.addAction(plg.menuentry, func(plg))
        return

    def save_plugin_data(self):
        """Save data for any plugins that need it"""

        data = {}
        for p in self.openplugins:
            plg = self.openplugins[p]
            #print (plg)
            d = plg.save_data()
            if d is not None:
                data[plg.name] = d
        #print (data)
        return data

    def clear_plugins(self):
        """Remove all open plugins"""

        keys = list(self.openplugins.keys())
        for p in keys:
            plg = self.openplugins[p]
            plg.main.deleteLater()
            del self.openplugins[plg.name]
            self.removeDockWidget(self.docks[plg.name])
            del self.docks[plg.name]
        return

    def preferences(self):
        """Preferences dialog"""

        opts = {}
        for k in core.defaults.keys():
            opts[k] = getattr(core,k)
        #opts['THEME'] = core.TEH
        dlg = widgets.PreferencesDialog(self, opts)
        dlg.exec_()
        return

    def online_documentation(self,event=None):
        """Open the online documentation"""

        link = 'https://github.com/dmnfarrell/snipgenie'
        link = 'https://snipgenie.readthedocs.io/'
        self.show_browser_tab(link, 'help')
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

    @Slot()
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
    finished = Signal()
    error = Signal(tuple)
    result = Signal(object)
    progress = Signal(str)

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
    """Class to provide a dialog for global options"""

    def __init__(self, parent=None):
        """Setup variables"""
        self.parent = parent
        self.kwds = {}
        genomes = []
        aligners = ['bwa','subread','minimap2']
        platforms = ['illumina','ont']
        separators = ['_','-','|',';','~','.']
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
