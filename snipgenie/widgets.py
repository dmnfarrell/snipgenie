# -*- coding: utf-8 -*-

"""
    Qt widgets for snpgenie.
    Created Jan 2020
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys, os, io
import numpy as np
import pandas as pd
import string
from .qt import *

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s
from . import tools, plotting

module_path = os.path.dirname(os.path.abspath(__file__))
iconpath = os.path.join(module_path, 'icons')

def dialogFromOptions(parent, opts, sections=None,
                      sticky='news', wrap=2, section_wrap=2):
    """Get Qt widgets dialog from a dictionary of options"""

    sizepolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    sizepolicy.setHorizontalStretch(1)
    sizepolicy.setVerticalStretch(0)

    style = '''
    QLabel {
        font-size: 12px;
    }
    QWidget {
        max-width: 130px;
        min-width: 30px;
        font-size: 14px;
    }
    QPlainTextEdit {
        max-height: 80px;
    }
    '''

    if sections == None:
        sections = {'options': opts.keys()}

    widgets = {}
    dialog = QWidget(parent)
    dialog.setSizePolicy(sizepolicy)

    l = QGridLayout(dialog)
    l.setSpacing(2)
    l.setAlignment(QtCore.Qt.AlignLeft)
    scol=1
    srow=1
    for s in sections:
        row=1
        col=1
        f = QGroupBox()
        f.setSizePolicy(sizepolicy)
        f.setTitle(s)
        #f.resize(50,100)
        #f.sizeHint()
        l.addWidget(f,srow,scol)
        gl = QGridLayout(f)
        gl.setAlignment(QtCore.Qt.AlignTop)
        srow+=1
        #gl.setSpacing(10)
        for o in sections[s]:
            label = o
            val = None
            opt = opts[o]
            if 'label' in opt:
                label = opt['label']
            val = opt['default']
            t = opt['type']
            lbl = QLabel(label)
            gl.addWidget(lbl,row,col)
            lbl.setStyleSheet(style)
            if t == 'combobox':
                w = QComboBox()
                w.addItems(opt['items'])
                #w.view().setMinListWidth(100)
                try:
                    w.setCurrentIndex(opt['items'].index(str(opt['default'])))
                except:
                    w.setCurrentIndex(0)
            elif t == 'entry':
                w = QLineEdit()
                w.setText(str(val))
            elif t == 'textarea':
                w = QPlainTextEdit()
                w.insertPlainText(str(val))
            elif t == 'slider':
                w = QSlider(QtCore.Qt.Horizontal)
                s,e = opt['range']
                w.setTickInterval(opt['interval'])
                w.setSingleStep(opt['interval'])
                w.setMinimum(s)
                w.setMaximum(e)
                w.setTickPosition(QSlider.TicksBelow)
                w.setValue(val)
            elif t == 'spinbox':
                if type(val) is float:
                    w = QDoubleSpinBox()
                else:
                    w = QSpinBox()
                w.setValue(val)
                if 'range' in opt:
                    min,max=opt['range']
                    w.setRange(min,max)
                    w.setMinimum(min)
                if 'interval' in opt:
                    w.setSingleStep(opt['interval'])
            elif t == 'checkbox':
                w = QCheckBox()
                w.setChecked(val)
            elif t == 'font':
                w = QFontComboBox()
                w.resize(w.sizeHint())
                w.setCurrentIndex(1)
            if 'width' in opt:
                h=20
                if 'height' in opt:
                    h=opt['height']
                w.setMinimumSize(opt['width'],h)
                w.resize(QtCore.QSize(opt['width'], h))

            #policy = dialog.sizePolicy()
            #policy.setVerticalStretch(1)
            #w.setSizePolicy(policy)

            col+=1
            gl.addWidget(w,row,col)
            w.setStyleSheet(style)
            widgets[o] = w
            #print (o, row, col)
            if col>=wrap:
                col=1
                row+=1
            else:
                col+=2
        if scol >= section_wrap:
            scol=1
        else:
            scol+=1
    return dialog, widgets

def getWidgetValues(widgets):
    """Get values back from a set of widgets"""

    kwds = {}
    for i in widgets:
        val = None
        if i in widgets:
            w = widgets[i]
            if type(w) is QLineEdit:
                try:
                    val = float(w.text())
                except:
                    val = w.text()
            elif type(w) is QPlainTextEdit:
                val = w.toPlainText()
            elif type(w) is QComboBox or type(w) is QFontComboBox:
                val = w.currentText()
            elif type(w) is QCheckBox:
                val = w.isChecked()
            elif type(w) is QSlider:
                val = w.value()
            elif type(w) in [QSpinBox,QDoubleSpinBox]:
                val = w.value()
            if val != None:
                kwds[i] = val
    kwds = kwds
    return kwds

def setWidgetValues(widgets, values):
    """Set values for a set of widgets from a dict"""

    kwds = {}
    for i in values:
        val = values[i]
        if i in widgets:
            #print (i, val, type(val))
            w = widgets[i]
            if type(w) is QLineEdit:
                w.setText(str(val))
            elif type(w) is QPlainTextEdit:
                w.insertPlainText(str(val))
            elif type(w) is QComboBox or type(w) is QFontComboBox:
                w.setCurrentIndex(1)
            elif type(w) is QCheckBox:
                w.setChecked(val)
            elif type(w) is QSlider:
                w.setValue(val)
            elif type(w) is QSpinBox:
                w.setValue(val)
    return

def addToolBarItems(toolbar, parent, items):
    """Populate toolbar from dict of items"""

    for i in items:
        if 'file' in items[i]:
            iconfile = os.path.join(iconpath,items[i]['file']+'.png')
            icon = QIcon(iconfile)
        else:
            icon = QIcon.fromTheme(items[i]['icon'])
        btn = QAction(icon, i, parent)
        btn.triggered.connect(items[i]['action'])
        if 'shortcut' in items[i]:
            btn.setShortcut(QKeySequence(items[i]['shortcut']))
        #btn.setCheckable(True)
        toolbar.addAction(btn)
    return toolbar

class MultipleInputDialog(QDialog):
    """Qdialog with multiple inputs"""
    def __init__(self, parent, options=None, title='Input', width=400, height=200):
        super(MultipleInputDialog, self).__init__(parent)
        self.values = None
        self.accepted = False
        self.setMinimumSize(width, height)
        self.setWindowTitle(title)
        dialog, self.widgets = dialogFromOptions(self, options)
        vbox = QVBoxLayout(self)
        vbox.addWidget(dialog)
        buttonbox = QDialogButtonBox(self)
        buttonbox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonbox.button(QDialogButtonBox.Ok).clicked.connect(self.accept)
        buttonbox.button(QDialogButtonBox.Cancel).clicked.connect(self.close)
        vbox.addWidget(buttonbox)
        self.show()
        return self.values

    def accept(self):
        self.values = getWidgetValues(self.widgets)
        self.accepted = True
        self.close()
        return

class ToolBar(QWidget):
    """Toolbar class"""
    def __init__(self, table, parent=None):
        super(ToolBar, self).__init__(parent)
        self.parent = parent
        self.table = table
        self.layout = QVBoxLayout()
        self.layout.setAlignment(QtCore.Qt.AlignTop)
        self.layout.setContentsMargins(2,2,2,2)
        self.setLayout(self.layout)
        self.createButtons()
        self.setMaximumWidth(40)
        return

    def createButtons(self):

        funcs = {'load':self.table.load, 'save':self.table.save,
                 'importexcel': self.table.load,
                 'copy':self.table.copy, 'paste':self.table.paste,
                 'plot':self.table.plot,
                 'transpose':self.table.pivot,
                 'pivot':self.table.pivot}
        icons = {'load': 'document-new', 'save': 'document-save-as',
                 'importexcel': 'x-office-spreadsheet',
                 'copy': 'edit-copy', 'paste': 'edit-paste',
                 'plot':'insert-image',
                 'transpose':'object-rotate-right',
                 'pivot': 'edit-undo',
                 }
        for name in funcs:
            self.addButton(name, funcs[name], icons[name])

    def addButton(self, name, function, icon):

        layout=self.layout
        button = QPushButton(name)
        button.setGeometry(QtCore.QRect(30,40,30,40))
        button.setText('')
        iconw = QIcon.fromTheme(icon)
        button.setIcon(QIcon(iconw))
        button.setIconSize(QtCore.QSize(20,20))
        button.clicked.connect(function)
        button.setMinimumWidth(30)
        layout.addWidget(button)

class BasicDialog(QDialog):
    """Qdialog for table operations interfaces"""
    def __init__(self, parent, table, title=None):

        super(BasicDialog, self).__init__(parent)
        self.parent = parent
        self.table = table
        self.df = table.model.df
        #self.app = self.parent.app
        self.setWindowTitle(title)
        self.createWidgets()
        self.setGeometry(QtCore.QRect(400, 300, 1000, 600))
        self.show()
        return

    def createWidgets(self):
        """Create widgets - override this"""

        cols = list(self.df.columns)

    def createButtons(self, parent):

        bw = self.button_widget = QWidget(parent)
        vbox = QVBoxLayout(bw)
        vbox.setAlignment(QtCore.Qt.AlignTop)
        button = QPushButton("Apply")
        button.clicked.connect(self.apply)
        vbox.addWidget(button)
        button = QPushButton("Update")
        button.clicked.connect(self.update)
        vbox.addWidget(button)
        button = QPushButton("Copy to clipboard")
        button.clicked.connect(self.copy_to_clipboard)
        vbox.addWidget(button)
        button = QPushButton("Close")
        button.clicked.connect(self.close)
        vbox.addWidget(button)
        return bw

    def apply(self):
        """Override this"""
        return

    def update(self):
        """Update the original table"""

        self.table.model.df = self.result.model.df
        self.table.refresh()
        self.close()
        return

    def copy_to_clipboard(self):
        """Copy result to clipboard"""

        df = self.result.model.df
        df.to_clipboard()
        return

    def close(self):
        self.destroy()
        return

class MergeDialog(BasicDialog):
    """Dialog to melt table"""
    def __init__(self, parent, table, df2, title='Merge Tables'):
        self.table = table
        self.df = table.model.df
        self.df2 = df2
        BasicDialog.__init__(self, parent, table, title)
        return

    def createWidgets(self):
        """Create widgets"""

        cols = self.df.columns
        cols2 = self.df2.columns
        ops = ['merge','concat']
        how = ['inner','outer','left','right']
        hbox = QHBoxLayout(self)
        main = QWidget(self)
        main.setMaximumWidth(300)
        hbox.addWidget(main)

        l = QVBoxLayout(main)
        w = self.ops_w = QComboBox(main)
        w.addItems(ops)
        l.addWidget(QLabel('Operation'))
        l.addWidget(w)
        w = self.lefton_w = QListWidget(main)
        w.setSelectionMode(QAbstractItemView.MultiSelection)
        w.addItems(cols)
        l.addWidget(QLabel('Left on'))
        l.addWidget(w)
        w = self.righton_w = QListWidget(main)
        w.setSelectionMode(QAbstractItemView.MultiSelection)
        w.addItems(cols2)
        l.addWidget(QLabel('Right on'))
        l.addWidget(w)

        w = self.leftindex_w = QCheckBox(main)
        w.setChecked(False)
        l.addWidget(QLabel('Use left index'))
        l.addWidget(w)
        w = self.rightindex_w = QCheckBox(main)
        w.setChecked(False)
        l.addWidget(QLabel('Use right index'))
        l.addWidget(w)

        w = self.how_w = QComboBox(main)
        w.addItems(how)
        l.addWidget(QLabel('How'))
        l.addWidget(w)

        w = self.left_suffw = QLineEdit('_1')
        l.addWidget(QLabel('Left suffix'))
        l.addWidget(w)
        w = self.right_suffw = QLineEdit('_2')
        l.addWidget(QLabel('Right suffix'))
        l.addWidget(w)

        from . import tables
        self.result = tables.DataFrameTable(self)
        hbox.addWidget(self.result)
        bf = self.createButtons(self)
        hbox.addWidget(bf)
        return

    def updateColumns(self):

        #self.df2 =
        cols2 = self.df2.columns
        return

    def apply(self):
        """Do the operation"""

        left_index = self.leftindex_w.isChecked()
        right_index = self.rightindex_w.isChecked()
        if left_index == True:
            lefton = None
        else:
            lefton = [i.text() for i in self.lefton_w.selectedItems()]
        if right_index == True:
            righton = None
        else:
            righton = [i.text() for i in self.righton_w.selectedItems()]
        how = self.how_w.currentText()
        op = self.ops_w.currentText()
        if op == 'merge':
            res = pd.merge(self.df, self.df2,
                            left_on=lefton,
                            right_on=righton,
                            left_index=left_index,
                            right_index=right_index,
                            how=how,
                            suffixes=(self.left_suffw .text(),self.right_suffw.text())
                            )
        else:
            res = pd.concat([self.df, self.df2])
        self.result.model.df = res
        self.result.refresh()
        return

class BaseOptions(object):
    """Class to generate widget dialog for dict of options"""
    def __init__(self, parent=None, opts={}, groups={}):
        """Setup variables"""

        self.parent = parent
        self.groups = groups
        self.opts = opts
        return

    def applyOptions(self):
        """Set the plot kwd arguments from the widgets"""

        self.kwds = getWidgetValues(self.widgets)
        return

    def apply(self):
        self.applyOptions()
        if self.callback != None:
            self.callback()
        return

    def showDialog(self, parent, wrap=2, section_wrap=2):
        """Auto create tk vars, widgets for corresponding options and
           and return the frame"""

        dialog, self.widgets = dialogFromOptions(parent, self.opts, self.groups,
                                wrap=wrap, section_wrap=section_wrap)
        return dialog

    def setWidgetValue(self, key, value):
        "Set a widget value"

        setWidgetValues(self.widgets, {key: value})
        self.applyOptions()
        return

    def increment(self, key, inc):
        """Increase the value of a widget"""

        new = self.kwds[key]+inc
        self.setWidgetValue(key, new)
        return

class DynamicDialog(QDialog):
    """Dynamic form using baseoptions"""

    def __init__(self, parent=None, options={}, groups=None, title='Dialog'):
        super(DynamicDialog, self).__init__(parent)
        self.setWindowTitle(title)
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.opts = BaseOptions(self, options, groups)
        dialog = self.opts.showDialog(self, wrap=1, section_wrap=1)
        layout.addWidget(dialog)
        buttonbox = QDialogButtonBox(self)
        buttonbox.addButton("Cancel", QDialogButtonBox.RejectRole)
        buttonbox.addButton("Ok", QDialogButtonBox.AcceptRole)
        self.connect(buttonbox, QtCore.SIGNAL("accepted()"), self, QtCore.SLOT("accept()"))
        self.connect(buttonbox, QtCore.SIGNAL("rejected()"), self, QtCore.SLOT("reject()"))
        layout.addWidget(buttonbox)
        return

    def get_values():
        """Get the widget values"""

        kwds = self.opts.kwds
        return kwds

class Editor(QTextEdit):
    def __init__(self, parent=None, **kwargs):
        super(Editor, self).__init__(parent, **kwargs)

    def zoom(self, delta):
        if delta < 0:
            self.zoomOut(1)
        elif delta > 0:
            self.zoomIn(1)

    def contextMenuEvent(self, event):

        menu = QMenu(self)
        copyAction = menu.addAction("Copy")
        clearAction = menu.addAction("Clear")
        zoominAction = menu.addAction("Zoom In")
        zoomoutAction = menu.addAction("Zoom Out")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == copyAction:
            self.copy()
        elif action == clearAction:
            self.clear()
        elif action == zoominAction:
            self.zoom(1)
        elif action == zoomoutAction:
            self.zoom(-1)

class PlainTextEditor(QPlainTextEdit):
    def __init__(self, parent=None, **kwargs):
        super(PlainTextEditor, self).__init__(parent, **kwargs)
        font = QFont("Monospace")
        font.setPointSize(10)
        font.setStyleHint(QFont.TypeWriter)
        self.setFont(font)
        return

    def zoom(self, delta):
        if delta < 0:
            self.zoomOut(1)
        elif delta > 0:
            self.zoomIn(1)

    def contextMenuEvent(self, event):

        menu = QMenu(self)
        copyAction = menu.addAction("Copy")
        clearAction = menu.addAction("Clear")
        zoominAction = menu.addAction("Zoom In")
        zoomoutAction = menu.addAction("Zoom Out")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == copyAction:
            self.copy()
        elif action == clearAction:
            self.clear()
        elif action == zoominAction:
            self.zoom(1)
        elif action == zoomoutAction:
            self.zoom(-1)

class TextViewer(QDialog):
    """Sequence records features viewer using dna_features_viewer"""
    def __init__(self, parent=None, text='', title='Text'):

        super(TextViewer, self).__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(QtCore.QRect(200, 200, 1000, 400))
        self.setMinimumHeight(150)
        self.add_widgets()
        self.ed.appendPlainText(text)
        return

    def add_widgets(self):
        """Add widgets"""

        l = QVBoxLayout(self)
        self.setLayout(l)
        font = QFont("Monospace")
        font.setPointSize(10)
        font.setStyleHint(QFont.TypeWriter)
        self.ed = ed = PlainTextEditor(self, readOnly=True)
        self.ed.setFont(font)
        l.addWidget(self.ed)
        self.show()

class FileViewer(QDialog):
    """Sequence records features viewer"""
    def __init__(self, parent=None, filename=None):

        #QDialog.__init__(self)
        super(FileViewer, self).__init__(parent)
        self.ed = ed = QPlainTextEdit(self, readOnly=True)
        #ed.setStyleSheet("font-family: monospace; font-size: 14px;")
        font = QFont("Monospace")
        font.setPointSize(10)
        font.setStyleHint(QFont.TypeWriter)
        self.ed.setFont(font)
        self.setWindowTitle('sequence features')
        self.setGeometry(QtCore.QRect(200, 200, 800, 800))
        #self.setCentralWidget(ed)
        l = QVBoxLayout(self)
        self.setLayout(l)
        l.addWidget(ed)
        self.show()

    def show_records(self, recs, format='genbank'):

        from Bio import SeqIO
        recs = SeqIO.to_dict(recs)
        if format == 'genbank':
            for r in recs:
                self.ed.appendPlainText(recs[r].format('genbank'))
        elif format == 'gff':
            tools.save_gff(recs,'temp.gff')
            f = open('temp.gff','r')
            for l in f.readlines():
                self.ed.appendPlainText(l)
        recnames = list(recs.keys())
        return

class PlotViewer(QDialog):
    """matplotlib plots widget"""
    def __init__(self, parent=None):

        super(PlotViewer, self).__init__(parent)
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
        self.setGeometry(QtCore.QRect(200, 200, 600, 600))
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        #self.show()
        #self.show_figure()
        return

    def show_figure(self, fig):

        from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
        import matplotlib.pyplot as plt

        #ax.plot(range(10))
        canvas = FigureCanvas(fig)
        self.grid.addWidget(canvas)
        self.fig = fig
        #self.ax = ax
        return

class BrowserViewer(QDialog):
    """matplotlib plots widget"""
    def __init__(self, parent=None):

        super(BrowserViewer, self).__init__(parent)
        self.add_widgets()
        return

    def add_widgets(self):
        """Add widgets"""

        layout = self.layout = QVBoxLayout()
        self.main = QWidget()
        vbox = QVBoxLayout(self.main)
        layout.addWidget(self.main)
        from PySide2.QtWebEngineWidgets import QWebEngineView
        self.browser = QWebEngineView()
        vbox = QVBoxLayout()
        self.setLayout(vbox)
        vbox.addWidget(self.browser)
        self.browser.setMinimumHeight(500)

        toolswidget = QWidget()
        toolswidget.setMaximumHeight(100)

        vbox.addWidget(toolswidget)
        l = QVBoxLayout(toolswidget)
        self.zoomslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(5)
        w.setMinimum(5)
        w.setMaximum(50)
        w.setValue(10)
        l.addWidget(w)
        w.valueChanged.connect(self.zoom)
        return

    def load_page(self, url):
        self.browser.setUrl(url)

    def zoom(self):
        zoom = self.zoomslider.value()/10
        self.browser.setZoomFactor(zoom)


class SimpleBamViewer(QDialog):
    """Sequence records features viewer using dna_features_viewer"""
    def __init__(self, parent=None, filename=None):

        super(SimpleBamViewer, self).__init__(parent)
        #self.setWindowTitle('Bam File View')
        self.setGeometry(QtCore.QRect(200, 200, 1000, 300))
        self.setMinimumHeight(150)
        self.fontsize = 8
        self.add_widgets()
        return

    def add_widgets(self):
        """Add widgets"""

        l = QVBoxLayout(self)
        self.setLayout(l)
        val=0
        navpanel = QWidget()
        navpanel.setMaximumHeight(60)
        l.addWidget(navpanel)
        bl = QHBoxLayout(navpanel)
        slider = QSlider(QtCore.Qt.Horizontal)
        slider.setTickPosition(slider.TicksBothSides)
        slider.setTickInterval(1000)
        slider.setPageStep(200)
        slider.setValue(1)
        #slider.sliderReleased.connect(self.value_changed)
        slider.valueChanged.connect(self.value_changed)
        self.slider = slider
        bl.addWidget(slider)

        backbtn = QPushButton('<')
        backbtn.setMaximumWidth(50)
        bl.addWidget(backbtn)
        backbtn.clicked.connect(self.prev_page)
        nextbtn = QPushButton('>')
        nextbtn.setMaximumWidth(50)
        bl.addWidget(nextbtn)
        nextbtn.clicked.connect(self.next_page)

        zoomoutbtn = QPushButton('-')
        zoomoutbtn.setMaximumWidth(50)
        bl.addWidget(zoomoutbtn)
        zoomoutbtn.clicked.connect(self.zoom_out)
        zoominbtn = QPushButton('+')
        zoominbtn.setMaximumWidth(50)
        bl.addWidget(zoominbtn)
        zoominbtn.clicked.connect(self.zoom_in)

        self.geneselect = QComboBox()
        self.geneselect.currentIndexChanged.connect(self.find_gene)
        bl.addWidget(self.geneselect)

        self.chromselect = QComboBox()
        self.chromselect.currentIndexChanged.connect(self.update_chrom)
        bl.addWidget(self.chromselect)
        self.textview = QTextEdit(readOnly=True)
        self.font = QFont("Monospace")
        self.font.setPointSize(self.fontsize)
        self.font.setStyleHint(QFont.TypeWriter)
        self.textview.setFont(self.font)
        #self.textview.setAutoFillBackground(True)
        #self.textview.setStyleSheet("QTextEdit { background-color: rgb(0, 0, 0); }")
        l.addWidget(self.textview)

        bottom = QWidget()
        bottom.setMaximumHeight(50)
        l.addWidget(bottom)
        bl2 = QHBoxLayout(bottom)
        self.loclbl = QLabel('')
        bl2.addWidget(self.loclbl)
        return

    def load_data(self, bam_file, ref_file, gb_file=None, vcf_file=None):
        """Load reference seq and get contig/chrom names"""

        self.ref_file = ref_file
        self.bam_file = bam_file
        self.gb_file = gb_file
        chromnames = plotting.get_fasta_names(ref_file)
        length = plotting.get_fasta_length(ref_file)
        if self.gb_file != None:
            df = tools.genbank_to_dataframe(gb_file)
            df.loc[df["gene"].isnull(),'gene'] = df.locus_tag
            genes = df.gene.unique()
            self.annot = df
            self.geneselect.addItems(genes)
            self.geneselect.setStyleSheet("QComboBox { combobox-popup: 0; }");
            self.geneselect.setMaxVisibleItems(10)
        self.chrom = chromnames[0]
        self.chromselect.addItems(chromnames)
        self.chromselect.setStyleSheet("QComboBox { combobox-popup: 0; }");
        self.chromselect.setMaxVisibleItems(10)
        sl = self.slider
        sl.setMinimum(1)
        sl.setMaximum(length)
        sl.setTickInterval(length/20)
        return

    def set_chrom(self, chrom):
        """Set the selected record which also updates the plot"""

        index = self.chromselect.findText(chrom)
        self.chromselect.setCurrentIndex(index)
        return

    def update_chrom(self, chrom=None):
        """Update after chromosome selection changed"""

        recname = self.chromselect.currentText()
        length = self.length = plotting.get_fasta_length(self.ref_file, key=chrom)
        sl = self.slider
        sl.setMinimum(1)
        sl.setMaximum(length)
        sl.setTickInterval(length/20)
        self.redraw()
        return

    def find_gene(self):
        """Go to selected gene if annotation present"""

        gene = self.geneselect.currentText()
        df = self.annot
        df = df[df.gene==gene]
        if len(df) == 0:
            return
        f = df.iloc[0]
        self.redraw(f.start)
        return

    def value_changed(self):
        """Callback for widgets"""

        length = self.length
        #r = self.view_range
        start = int(self.slider.value())
        self.redraw(start)
        return

    def redraw(self, xstart=1):
        """Plot the features"""

        h = 5
        length = self.length
        if xstart<0:
            xstart=1
        ref = self.ref_file
        txt = tools.samtools_tview(self.bam_file, self.chrom, xstart, 500, ref, 'H')
        self.textview.setText(txt)
        self.loclbl.setText(str(xstart))
        return

    def prev_page(self):
        """ """
        start = int(self.slider.value())
        start -= 500
        if start<0:
            return
        self.redraw(start)
        self.slider.setValue(start)
        return

    def next_page(self):
        """ """
        start = int(self.slider.value())
        start += 500
        self.redraw(start)
        self.slider.setValue(start)
        return

    def zoom_in(self):
        """Zoom in"""

        self.fontsize = self.fontsize+1
        self.font.setPointSize(self.fontsize)
        self.textview.setFont(self.font)
        self.redraw()
        return

    def zoom_out(self):
        """Zoom out"""

        if self.fontsize <= 2:
            return
        self.fontsize = self.fontsize-1
        self.font.setPointSize(self.fontsize)
        self.textview.setFont(self.font)
        self.redraw()
        return

class GraphicalBamViewer(QDialog):
    """Alignment viewer with pylab"""
    def __init__(self, parent=None, filename=None):

        super(BamViewer, self).__init__(parent)
        #self.setWindowTitle('Bam File View')
        self.setGeometry(QtCore.QRect(200, 200, 1000, 300))
        self.setMinimumHeight(150)
        self.add_widgets()
        return

    def add_widgets(self):
        """Add widgets"""

        l = QVBoxLayout(self)
        self.setLayout(l)
        val=0
        navpanel = QWidget()
        navpanel.setMaximumHeight(60)
        l.addWidget(navpanel)
        bl = QHBoxLayout(navpanel)
        slider = QSlider(QtCore.Qt.Horizontal)
        slider.setTickPosition(slider.TicksBothSides)
        slider.setTickInterval(1000)
        slider.setPageStep(200)
        slider.setValue(1)
        #slider.sliderReleased.connect(self.value_changed)
        slider.valueChanged.connect(self.value_changed)
        self.slider = slider
        bl.addWidget(slider)

        zoomoutbtn = QPushButton('-')
        zoomoutbtn.setMaximumWidth(50)
        bl.addWidget(zoomoutbtn)
        zoomoutbtn.clicked.connect(self.zoom_out)
        zoominbtn = QPushButton('+')
        zoominbtn.setMaximumWidth(50)
        bl.addWidget(zoominbtn)
        zoominbtn.clicked.connect(self.zoom_in)

        self.chromselect = QComboBox()
        self.chromselect.currentIndexChanged.connect(self.update_chrom)
        bl.addWidget(self.chromselect)

        #add plot axes
        from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(15,2))
        spec = fig.add_gridspec(ncols=1, nrows=4)
        self.ax1 = fig.add_subplot(spec[0, 0])
        self.ax2 = fig.add_subplot(spec[1:3, 0])
        self.ax3 = fig.add_subplot(spec[3, 0])
        self.canvas = FigureCanvas(fig)
        self.fig = fig
        l.addWidget(self.canvas)

        bottom = QWidget()
        bottom.setMaximumHeight(50)
        l.addWidget(bottom)
        bl2 = QHBoxLayout(bottom)
        self.loclbl = QLabel('')
        bl2.addWidget(self.loclbl)
        return

    def load_data(self, bam_file, ref_file, gb_file=None, vcf_file=None):
        """Load reference seq and get contig/chrom names"""

        self.ref_file = ref_file
        self.bam_file = bam_file
        self.gb_file = gb_file
        chromnames = plotting.get_fasta_names(ref_file)
        length = plotting.get_fasta_length(ref_file)
        self.chrom = chromnames[0]
        self.chromselect.addItems(chromnames)
        self.chromselect.setStyleSheet("QComboBox { combobox-popup: 0; }");
        self.chromselect.setMaxVisibleItems(10)
        sl = self.slider
        sl.setMinimum(1)
        sl.setMaximum(length)
        sl.setTickInterval(length/20)
        return

    def set_chrom(self, chrom):
        """Set the selected record which also updates the plot"""

        index = self.chromselect.findText(chrom)
        self.chromselect.setCurrentIndex(index)
        return

    def update_chrom(self, chrom=None):
        """Update after chromosome selection changed"""

        recname = self.chromselect.currentText()
        length = self.length = plotting.get_fasta_length(self.ref_file, key=chrom)
        sl = self.slider
        sl.setMinimum(1)
        sl.setMaximum(length)
        sl.setTickInterval(length/20)
        self.redraw()
        return

    def value_changed(self):
        """Callback for widgets"""

        length = self.length
        r = self.view_range
        start = int(self.slider.value())
        end = int(start+r)
        if end > length:
            return
        self.redraw(start, end)
        return

    def zoom_in(self):
        """Zoom in"""

        #length = len(self.rec.seq)
        fac = 1.2
        r = int(self.view_range/fac)
        start = int(self.slider.value())
        end = start + r
        #if end > length:
        #    end=length
        self.redraw(start, end)
        return

    def zoom_out(self):
        """Zoom out"""

        #length = len(self.rec.seq)
        fac = 1.2
        r = int(self.view_range*fac)
        start = int(self.slider.value())
        end = start + r
        #if end > length:
        #    end=length
        #    start = start-r
        self.redraw(start, end)
        return

    def redraw(self, xstart=1, xend=2000):
        """Plot the features"""

        h = 5
        self.ax1.clear()
        self.ax2.clear()
        length = self.length
        if xstart<0:
            xstart=1
        if xend <= 0:
            xend = xstart+2000
        if xend-xstart > 10000:
            xend = xstart+10000
        if xend > length:
            xend = length

        #print (start, end)
        cov = plotting.get_coverage(self.bam_file, self.chrom, xstart, xend)
        plotting.plot_coverage(cov,ax=self.ax1,xaxis=False)
        plotting.plot_bam_alignment(self.bam_file, self.chrom, xstart, xend, ax=self.ax2)
        if self.gb_file != None:
            recs = tools.gb_to_records(self.gb_file)
            plotting.plot_features(recs[0], self.ax3, xstart=xstart, xend=xend)

        self.canvas.draw()
        self.view_range = xend-xstart
        self.loclbl.setText(str(xstart)+'-'+str(xend))
        return

class TreeViewer(QDialog):
    """Phylogeny viewer with toytree"""
    def __init__(self, parent=None, filename=None):

        super(TreeViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 1000, 300))
        self.setMinimumHeight(150)
        self.add_widgets()
        self.create_menu(self)
        self.tree = None
        self.width = 400
        self.height = 500
        self.ts = ''
        import toytree
        self.colors = {}
        self.default_style = {
            "layout":'r',
            "edge_type": 'p',
            "edge_style": {
                "stroke": 'black',
                "stroke-width": 2,
            },
            "tip_labels": True,
            "tip_labels_align": True,
            "tip_labels_colors": 'black',
            "tip_labels_style": {
                "font-size": "14px"
            },
            "node_labels": False,
            "node_sizes": 10,
            "node_colors": toytree.colors[2],
            "node_markers":"c",
            "use_edge_lengths":True,
        }
        self.style = self.default_style.copy()

        self.test_tree(10)
        return

    def test_tree(self, n=None):
        """Load a test tree"""

        import toytree
        if n==None:
            n, ok = QInputDialog().getInt(self, "Test tree",
                                                 "Nodes:", 10)
            if not ok:
                return
        self.set_tree(self.random_tree(n=n))
        self.height = 200+self.tree.ntips*10
        self.update()
        return

    def random_tree(self, n=12):
        """Make a random tree"""

        import toytree
        tre = toytree.rtree.coaltree(n)
        ## assign random edge lengths and supports to each node
        for node in tre.treenode.traverse():
            node.dist = np.random.exponential(1)
            node.support = int(np.random.uniform(50, 100))
        return tre

    def saveData(self):
        """Save layers"""

        data = tools.get_attributes(self)
        data['tree'] = self.tree
        return data

    def loadData(self, data):
        """Load saved layers"""

        try:
            self.set_tree(data['tree'])
            tools.set_attributes(self, data)
        except:
            pass
        self.update()
        return

    def create_menu(self, parent):
        """Menu bar"""

        self.menubar = QMenuBar(parent)
        self.file_menu = QMenu('File', parent)
        self.file_menu.addAction('Import Tree', self.load_tree)
        self.file_menu.addAction('Import MultiTree', self.load_multitree)
        self.file_menu.addAction('Load Test Tree', self.test_tree)
        self.file_menu.addAction('Show Newick', self.show_newick)
        self.file_menu.addAction('Export Image', self.export_image)
        self.menubar.addMenu(self.file_menu)
        self.tree_menu = QMenu('Tree', parent)
        self.tree_menu.addAction('Show Unrooted', self.unroot_tree)
        self.tree_menu.addAction('Reset Format', self.reset_style)
        self.menubar.addMenu(self.tree_menu)

        return

    def add_widgets(self):
        """Add widgets"""

        layout = self.layout = QVBoxLayout()
        self.main = QWidget()
        vbox = QVBoxLayout(self.main)
        layout.addWidget(self.main)
        from PySide2.QtWebEngineWidgets import QWebEngineView
        self.browser = QWebEngineView()
        vbox = QVBoxLayout()
        self.setLayout(vbox)
        vbox.addWidget(self.browser)
        self.browser.setMinimumHeight(500)

        toolswidget = QWidget()
        toolswidget.setMaximumHeight(200)

        vbox.addWidget(toolswidget)
        l = QVBoxLayout(toolswidget)
        self.zoomslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(5)
        w.setMinimum(5)
        w.setMaximum(50)
        w.setValue(10)
        l.addWidget(w)
        w.valueChanged.connect(self.zoom)
        btn = QPushButton('Set Format')
        l.addWidget(btn)
        btn.clicked.connect(self.tree_style_options)
        t = self.tipitems = QTreeWidget()
        t.setHeaderItem(QTreeWidgetItem(["name","visible"]))
        t.setColumnWidth(0, 200)
        t.setSelectionMode(QAbstractItemView.ExtendedSelection)
        t.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        t.customContextMenuRequested.connect(self.show_tree_menu)
        l.addWidget(t)

        return

    def show_tree_menu(self, pos):
        """Show right cick tree menu"""

        item = self.tipitems.itemAt( pos )
        menu = QMenu(self.tipitems)
        colorAction = menu.addAction("Set Color")
        nodecolorAction = menu.addAction("Set Node Color")
        rootAction = menu.addAction("Root On")
        dropAction = menu.addAction("Drop Tips")
        action = menu.exec_(self.tipitems.mapToGlobal(pos))
        if action == rootAction:
            self.root_tree()
        elif action == colorAction:
            self.set_color()
        elif action == nodecolorAction:
            self.set_color('node')
        elif action == dropAction:
            self.drop_tips()

    def load_tree(self, filename=None):
        import toytree
        if filename == None:
            options = QFileDialog.Options()
            filter = "newick files (*.newick);;All files (*.*)"
            filename, _ = QFileDialog.getOpenFileName(self,"Open tree file",
                                        "",filter=filter,selectedFilter =filter, options=options)
            if not filename:
                return
        self.set_tree(toytree.tree(filename))
        return

    def set_tree(self, tree):
        """Set a new tree"""

        self.tree = tree
        self.colors = {}
        self.style['tip_labels_colors'] = 'black'
        self.tipitems.clear()
        for t in self.tree.get_tip_labels():
            item = QTreeWidgetItem(self.tipitems)
            item.setCheckState(1, QtCore.Qt.Checked)
            item.setText(0, t)

        return

    def load_multitree(self):

        import toytree
        options = QFileDialog.Options()
        filter = "newick files (*.newick);;All files (*.*)"
        filename, _ = QFileDialog.getOpenFileName(self,"Open tree file",
                                    "",filter=filter,selectedFilter =filter, options=options)
        if not filename:
            return
        self.tree = toytree.mtree(filename)
        self.update()
        return

    def update(self):
        """Update the plot"""

        import toytree
        import toyplot

        if self.tree==None:
            return
        #set colors
        print (type(self.tree))
        if type(self.tree) is toytree.Multitree.MultiTree:
            self.update_multitree()
            return
        colorlist = [self.colors[tip] if tip in self.colors else "black" for tip in self.tree.get_tip_labels()]
        self.style['tip_labels_colors'] = colorlist
        if self.ts != '':
            style = {}
        else:
            style = self.style
        canvas,axes,mark = self.tree.draw(ts=self.ts,
                        width=self.width,
                        height=self.height,
                        scalebar=True, **style)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            self.browser.setHtml(html)
        self.canvas = canvas
        return

    def update_multitree(self):
        import toyplot
        style = self.style
        canvas,axes,mark = self.tree.draw(ncols=3,nrows=3, ts=self.ts,
                        width=self.width,
                        height=self.height*3, **style)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            self.browser.setHtml(html)
        self.canvas = canvas
        return

    def show_newick(self):

        txt = self.tree.newick
        ed = TextViewer(self,txt)
        return

    def root_tree(self):

        item = self.tipitems.selectedItems()[0]
        row = self.tipitems.selectedIndexes()[0].row()
        name = item.text(0)
        self.tree = self.tree.root(name).ladderize()
        self.update()
        return

    def unroot_tree(self):
        self.tree = self.tree.unroot().ladderize()
        self.update()
        return

    def export_image(self):
        """Save tree as image"""

        options = QFileDialog.Options()
        filter = "png files (*.png);;pdf files (*.pdf);;All files (*.*)"
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                    "",filter=filter,selectedFilter =filter, options=options)
        if not filename:
            return

        ext = os.path.splitext(filename)
        print (ext)
        from toyplot import png
        png.render(self.canvas, filename, width=(4, "inches"))
        return

    def zoom(self):
        zoom = self.zoomslider.value()/10
        self.browser.setZoomFactor(zoom)

    def clustermap(self):
        #canvas = toyplot.Canvas(width=600, height=800)
        #axes = canvas.cartesian()
        #y = np.linspace(0, 1, 20) ** 2
        #axes.plot(y)
        #matrix = np.random.normal(loc=1.0, size=(20, 20))
        #canvas.matrix(matrix, label="A matrix")
        return

    def tree_style_options(self):

        fonts = ['%spx' %i for i in range (6,28)]
        markers = ['o','s','d','^','>','oo']
        nlabels = ['','idx','support']
        tip_labels_style = self.style['tip_labels_style']
        tree_styles = [None,'n','s','c','o','p','d']
        opts = {'ts': {'type':'combobox','default':self.ts,'items':tree_styles},
                'layout': {'type':'combobox','default':self.style['layout'],'items':['r','d','c']},
                'edge_type': {'type':'combobox','default':self.style['edge_type'],'items':['p','b','c']},
                'use_edge_lengths':{'type':'checkbox','default':self.style['use_edge_lengths'] },
                'tip_labels':{'type':'checkbox','default':self.style['tip_labels'] },
                'tip_labels_align':{'type':'checkbox','default':self.style['tip_labels_align'] },
                'node_labels':{'type':'combobox','default':self.style['node_labels'],'items': nlabels},
                'node_sizes':{'type':'spinbox','default':self.style['node_sizes'],'range':(2,20),'interval':1},
                'node_markers': {'type':'combobox','default':self.style['node_markers'],'items':markers},
                'font_size':{'type':'combobox','default':tip_labels_style['font-size'],'items':fonts},
                'width':{'type':'entry','default':self.width},
                'height':{'type':'entry','default':self.height,},
                }

        dlg = MultipleInputDialog(self, opts, title='Tree Style', width=300)
        dlg.exec_()
        if not dlg.accepted:
            return False
        kwds = dlg.values
        self.set_style(kwds)
        self.update()
        return

    def set_style(self, kwds):

        omit=['width','height','font_size','ts']
        for k in kwds:
            if k not in omit:
                self.style[k] = kwds[k]
        if kwds['node_labels'] == '':
            self.style['node_labels'] = False
        self.style['tip_labels_style']['font-size'] = kwds['font_size']
        self.ts = kwds['ts']
        self.width = kwds['width']
        self.height = kwds['height']
        return

    def reset_style(self):

        self.style = self.default_style
        self.colors = {}
        print (self.style)
        self.update()

    def set_color(self, kind='text'):

        items = self.tipitems.selectedItems()
        names = [i.text(0) for i in items]
        qcolor = QColorDialog.getColor()
        for item in items:
            item.setBackground(0 , qcolor)
        for name in names:
            if kind == 'text':
                self.colors[name] = qcolor.name()
            elif kind == 'node':
                self.node_colors[name] = qcolor.name()
        self.update()
        return

    def drop_tips(self):

        items = self.tipitems.selectedItems()
        names = [i.text(0) for i in items]
        #for name in names:
        self.tree = self.tree.drop_tips(names=names).ladderize()
        self.update()
        return
