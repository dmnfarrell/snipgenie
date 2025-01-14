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

import sys, os, io, platform, time
import pickle
import numpy as np
import pandas as pd
import string
import pylab as plt
import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from .qt import *
from . import tables

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s
from . import core, tools, plotting

module_path = os.path.dirname(os.path.abspath(__file__))
iconpath = os.path.join(module_path, 'icons')

def create_button(parent, name, function, iconname=None, iconsize=20,
                 tooltip=None):
    """Create a button for a function and optional icon.
        Returns:
            the button widget
    """

    button = QPushButton(parent)
    #button.setGeometry(QtCore.QRect(40,40,40,40))
    button.setText(name)
    iconfile = os.path.join(iconpath,iconname)
    icon = QIcon(iconfile)
    button.setIcon(QIcon(icon))
    button.setIconSize(QtCore.QSize(iconsize,iconsize))
    button.clicked.connect(function)
    #button.setMinimumWidth(20)
    if tooltip != None:
        button.setToolTip(tooltip)
    return button

def dialogFromOptions(parent, opts, sections=None,
                      wrap=2, section_wrap=4,
                      style=None):
    """
    Get Qt widgets dialog from a dictionary of options.
    Args:
        opts: options dictionary
        sections:
        section_wrap: how many sections in one row
        style: stylesheet css if required
    """

    sizepolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    sizepolicy.setHorizontalStretch(0)
    sizepolicy.setVerticalStretch(0)

    if style == None:
        style = '''
        QWidget {
            font-size: 12px;
        }
        QPlainTextEdit {
            max-height: 80px;
        }
        QComboBox {
            combobox-popup: 0;
            max-height: 30px;
            max-width: 150px;
        }
        '''

    if sections == None:
        sections = {'options': opts.keys()}

    widgets = {}
    dialog = QWidget(parent)
    dialog.setSizePolicy(sizepolicy)

    l = QGridLayout(dialog)
    l.setSpacing(1)
    l.setAlignment(QtCore.Qt.AlignTop)
    scol=1
    srow=1
    for s in sections:
        row=srow
        col=1
        f = QWidget()
        f.resize(50,100)
        f.sizeHint()
        l.addWidget(f,row,scol)
        gl = QGridLayout(f)
        gl.setAlignment(QtCore.Qt.AlignTop)
        gl.setSpacing(5)
        for o in sections[s]:
            label = o
            val = None
            opt = opts[o]
            if 'label' in opt:
                label = opt['label']
            val = opt['default']
            t = opt['type']
            lbl = QLabel(label)
            lbl.setMinimumWidth(150)
            gl.addWidget(lbl,row,col)
            lbl.setStyleSheet(style)
            if t == 'combobox':
                w = QComboBox()
                w.addItems(opt['items'])
                index = w.findText(val)
                if index != -1:
                    w.setCurrentIndex(index)
                if 'editable' in opt:
                     w.setEditable(True)
                if 'width' in opt:
                    w.setMinimumWidth(opt['width'])
                    w.resize(opt['width'], 20)
                w.view().setMinimumWidth(120)
                w.setMaxVisibleItems(12)
            elif t == 'list':
                w = QListWidget()
                w.setSelectionMode(QAbstractItemView.MultiSelection)
                w.addItems(opt['items'])
            elif t == 'entry':
                w = QLineEdit()
                w.setText(str(val))
                if 'width' in opt:
                    w.setMaximumWidth(opt['width'])
                    w.resize(opt['width'], 20)
            elif t == 'textarea':
                w = QPlainTextEdit()
                #w.setSizePolicy(sizepolicy)
                w.insertPlainText(str(val))
                w.setMaximumHeight(100)
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
                w = QSpinBox()
                w.setValue(val)
                if 'range' in opt:
                    min, max = opt['range']
                    min = int(min)
                    max = int(max)
                    w.setRange(min,max)
                    w.setMaximum(max)
                    w.setMinimum(min)
                if 'interval' in opt:
                    w.setSingleStep(opt['interval'])
            elif t == 'doublespinbox':
                w = QDoubleSpinBox()
                w.setValue(val)
                if 'range' in opt:
                    min, max = opt['range']
                    w.setRange(min,max)
                    w.setMinimum(min)
                if 'interval' in opt:
                    w.setSingleStep(opt['interval'])
            elif t == 'checkbox':
                w = QCheckBox()
                w.setChecked(val)
            elif t == 'font':
                w = QFontComboBox()
                index = w.findText(val)
                #w.resize(w.sizeHint())
                w.setCurrentIndex(index)
            elif t == 'dial':
                w = QDial()
                if 'range' in opt:
                    min, max = opt['range']
                    w.setMinimum(min)
                    w.setMaximum(max)
                w.setValue(val)
            elif t == 'colorbutton':
                w = ColorButton()
                w.setColor(val)
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
            srow+=2
        else:
            scol+=1
    return dialog, widgets

def getWidgetValue(w):
    """Get value from any kind of widget"""

    val = None
    if type(w) is QLineEdit:
        val = w.text()
    elif type(w) is QPlainTextEdit:
        val = w.toPlainText()
    elif type(w) is QComboBox or type(w) is QFontComboBox:
        val = w.currentText()
    elif type(w) is QListWidget:
        val = [i.text() for i in w.selectedItems()]
    elif type(w) is QCheckBox:
        val = w.isChecked()
    elif type(w) is QSlider:
        val = w.value()
    elif type(w) in [QSpinBox,QDoubleSpinBox]:
        val = w.value()
    elif type(w) is ColorButton:
        val = w.color()
    return val

def getWidgetValues(widgets):
    """Get values back from a set of widgets"""

    kwds = {}
    for i in widgets:
        val = None
        if i in widgets:
            w = widgets[i]
            val = getWidgetValue(w)
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
                index = w.findText(val)
                w.setCurrentIndex(index)
            elif type(w) is QCheckBox:
                w.setChecked(val)
            elif type(w) is QSlider:
                w.setValue(val)
            elif type(w) is QSpinBox:
                w.setValue(int(val))
            elif type(w) is QDoubleSpinBox:
                w.setValue(val)
    return

def addToolBarItems(toolbar, parent, items):
    """Populate toolbar from dict of items"""

    for i in items:
        if 'file' in items[i]:
            iconfile = os.path.join(iconpath,items[i]['file'])
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

def get_fonts():
     """Get the current list of system fonts"""

     import matplotlib.font_manager
     l = matplotlib.font_manager.findSystemFonts(fontext='ttf')
     fonts = []
     for fname in l:
        try: fonts.append(matplotlib.font_manager.FontProperties(fname=fname).get_name())
        except RuntimeError: pass
     fonts = list(set(fonts))
     fonts.sort()
     return fonts

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

class ColorButton(QPushButton):
    '''
    Custom Qt Widget to show a chosen color.

    Left-clicking the button shows the color-chooser, while
    right-clicking resets the color to None (no-color).
    '''

    colorChanged = Signal(object)

    def __init__(self, *args, color=None, **kwargs):
        super(ColorButton, self).__init__(*args, **kwargs)

        self._color = None
        self._default = color
        self.pressed.connect(self.onColorPicker)

        # Set the initial/default state.
        self.setColor(self._default)

    def setColor(self, color):

        if color != self._color:
            self._color = color
            self.colorChanged.emit(color)

        if self._color:
            self.setStyleSheet("background-color: %s;" % self._color)
        else:
            self.setStyleSheet("")
        return

    def color(self):
        return self._color

    def onColorPicker(self):
        '''
        Show color-picker dialog to select color.
        Qt will use the native dialog by default.
        '''
        dlg = QColorDialog(self)
        if self._color:
            dlg.setCurrentColor(QColor(self._color))

        if dlg.exec_():
            self.setColor(dlg.currentColor().name())

    def mousePressEvent(self, e):
        if e.button() == QtCore.Qt.RightButton:
            self.setColor(self._default)

        return super(ColorButton, self).mousePressEvent(e)

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
        how = ['left','inner','right','outer']
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

class PreferencesDialog(QDialog):
    """Preferences dialog from config parser options"""

    def __init__(self, parent, options={}):

        super(PreferencesDialog, self).__init__(parent)
        self.parent = parent
        self.setWindowTitle('Preferences')
        self.resize(400, 200)
        self.setGeometry(QtCore.QRect(300,300, 600, 200))
        self.setMaximumWidth(600)
        self.setMaximumHeight(300)
        self.createWidgets(options)
        self.show()
        return

    def createWidgets(self, options):
        """create widgets"""

        import pylab as plt

        colormaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        timeformats = ['%m/%d/%Y','%d/%m/%Y','%d/%m/%y',
                '%Y/%m/%d','%y/%m/%d','%Y/%d/%m',
                '%d-%b-%Y','%b-%d-%Y',
                '%Y-%m-%d %H:%M:%S','%Y-%m-%d %H:%M',
                '%d-%m-%Y %H:%M:%S','%d-%m-%Y %H:%M',
                '%Y','%m','%d','%b']
        plotstyles = [''] + plt.style.available
        themes = QStyleFactory.keys() + ['dark','light']

        self.opts = {
                'FONT':{'type':'font','default':options['FONT'],'label':'Font'},
                'FONTSIZE':{'type':'spinbox','default':options['FONTSIZE'],'range':(5,40),
                            'interval':1,'label':'Font Size'},
                'TIMEFORMAT':{'type':'combobox','default':options['TIMEFORMAT'],
                            'items':timeformats,'label':'Date/Time format'},
                'PLOTSTYLE':{'type':'combobox','default':options['PLOTSTYLE'],
                            'items':plotstyles,'label':'Plot Style'},
                'DPI':{'type':'entry','default':options['DPI'],#'range':(20,300),'interval':10,
                        'label':'Plot DPI'},
                'ICONSIZE':{'type':'spinbox','default':options['ICONSIZE'],'range':(16,64), 'label':'Icon Size'},
                }
        sections = {'formats':['FONT','FONTSIZE','TIMEFORMAT','ICONSIZE','PLOTSTYLE','DPI']
                    }

        dialog, self.widgets = dialogFromOptions(self, self.opts, sections)

        self.layout = QVBoxLayout(self)
        self.layout.addWidget(dialog)
        dialog.setFocus()
        bw = self.createButtons(self)
        self.layout.addWidget(bw)
        return

    def createButtons(self, parent):

        bw = self.button_widget = QWidget(parent)
        vbox = QHBoxLayout(bw)
        button = QPushButton("Apply")
        button.clicked.connect(self.apply)
        vbox.addWidget(button)
        button = QPushButton("Reset")
        button.clicked.connect(self.reset)
        vbox.addWidget(button)
        button = QPushButton("Close")
        button.clicked.connect(self.close)
        vbox.addWidget(button)
        return bw

    def apply(self):
        """Apply options to current table"""

        kwds = getWidgetValues(self.widgets)
        core.FONT = kwds['FONT']
        core.FONTSIZE = kwds['FONTSIZE']
        core.TIMEFORMAT = kwds['TIMEFORMAT']
        core.PLOTSTYLE = kwds['PLOTSTYLE']
        core.DPI = kwds['DPI']
        core.ICONSIZE = kwds['ICONSIZE']

        self.parent.refresh()
        self.parent.apply_settings()
        return

    def updateWidgets(self, kwds=None):
        """Update widgets from stored or supplied kwds"""

        if kwds == None:
            kwds = self.kwds
        for k in kwds:
            setWidgetValues(self.widgets, {k: kwds[k]})
        return

    def setDefaults(self):
        """Populate default kwds dict"""

        self.kwds = {}
        for o in self.opts:
            self.kwds[o] = core.defaults[o]
        return

    def reset(self):
        """Reset to defaults"""

        self.setDefaults()
        self.updateWidgets()
        self.apply()
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

    def showDialog(self, parent, wrap=2, section_wrap=2, style=None):
        """Auto create tk vars, widgets for corresponding options and
           and return the frame"""

        dialog, self.widgets = dialogFromOptions(parent, self.opts, self.groups,
                                wrap=wrap, section_wrap=section_wrap, style=style)
        return dialog

    def setWidgetValue(self, key, value):
        "Set a widget value"

        setWidgetValues(self.widgets, {key: value})
        self.applyOptions()
        return

    def updateWidgets(self, kwds):

        for k in kwds:
            setWidgetValues(self.widgets, {k: kwds[k]})
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
    def __init__(self, parent=None, fontsize=12, **kwargs):
        super(Editor, self).__init__(parent, **kwargs)
        font = QFont("Monospace")
        font.setPointSize(fontsize)
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

    def insert(self, txt):

        self.insertPlainText(txt)
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())
        return

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
    """Plain text viewer"""
    def __init__(self, parent=None, text='', width=200, height=400, title='Text'):

        super(TextViewer, self).__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(QtCore.QRect(200, 200, width, height))
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

class PlotWidget(FigureCanvas):
    """Basic mpl plot view."""

    def __init__(self, parent=None, figure=None, dpi=100, hold=False):

        if figure == None:
            figure = Figure()
        super(PlotWidget, self).__init__(figure)
        self.setParent(parent)
        if figure is not None:
            self.fig = figure
        else:
            self.fig = Figure(dpi=dpi)
        self.canvas = FigureCanvas(self.figure)
        #self.ax = self.figure.add_subplot(111)

    def clear(self):
        """Clear plot"""

        self.fig.clear()
        self.canvas.draw()
        return

    def set_figure(self, fig):
        """Set the figure if we have plotted elsewhere"""

        self.clear()
        self.fig = fig
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        return

class BasePlotViewer(QWidget):
    """matplotlib plots widget base class"""
    def __init__(self, parent=None, title=''):

        super(BasePlotViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 900, 600))
        self.main = QSplitter()
        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.main)
        self.create_figure()
        #self.create_controls()
        self.setWindowTitle(title)
        return

    def create_figure(self, fig=None):
        """Create canvas and figure"""

        import matplotlib.pyplot as plt
        #ax.plot(range(10))
        if fig == None:
            fig, ax = plt.subplots(1,1, figsize=(7,5), dpi=120)
            self.ax = ax
        if hasattr(self, 'canvas'):
            self.layout().removeWidget(self.canvas)
        canvas = FigureCanvas(fig)
        left = QWidget()
        self.main.addWidget(left)
        l = QVBoxLayout()
        left.setLayout(l)
        l.addWidget(canvas)
        self.toolbar = NavigationToolbar(canvas, self)
        l.addWidget(self.toolbar)
        self.fig = fig
        self.canvas = canvas
        iconfile = os.path.join(iconpath,'reduce')
        a = QAction(QIcon(iconfile), "Reduce elements",  self)
        a.triggered.connect(lambda: self.zoom(zoomin=False))
        self.toolbar.addAction(a)
        iconfile = os.path.join(iconpath,'enlarge')
        a = QAction(QIcon(iconfile), "Enlarge elements",  self)
        a.triggered.connect(lambda: self.zoom(zoomin=True))
        self.toolbar.addAction(a)
        return

    def set_figure(self, fig):
        """Set the figure if we have plotted elsewhere"""

        self.clear()
        self.create_figure(fig)
        #self.ax = fig.ax
        self.canvas.draw()
        return

    def clear(self):
        """Clear plot"""

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.canvas.draw()
        return

class PlotOptions(BaseOptions):
    """Class to provide a dialog for plot options"""

    def __init__(self, parent=None):
        """Setup variables"""

        self.parent = parent
        self.kwds = {}
        kinds = ['','bar','barh','hist','scatter','line','heatmap','pie','box']
        scales = ['linear','log']
        style_list = plt.style.available
        markers = ['','o','.','^','v','>','<','s','+','x','p','d','h','*']
        linestyles = ['-','--','-.',':']
        fonts = get_fonts()
        if 'Windows' in platform.platform():
            defaultfont = 'Arial'
        else:
            defaultfont = 'FreeSans'
        colormaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        self.groups = {'general':['kind','subplots','grid','sharex','bins','linewidth','linestyle',
                       'marker','ms','alpha','colormap'],
                       'format' :['title','xlabel','ylabel','showxlabels','showylabels',
                       'style','font','fontsize']
                       }
        self.opts = {
                    'kind':{'type':'combobox','default':'','items':kinds},
                    'grid':{'type':'checkbox','default':0,'label':'show grid'},
                    'subplots':{'type':'checkbox','default':0,'label':'subplots'},
                    'sharex':{'type':'checkbox','default':0,'label':'share x'},
                    'bins':{'type':'spinbox','default':20,'width':5},
                    'marker':{'type':'combobox','default':'o','items': markers},
                    'linestyle':{'type':'combobox','default':'-','items': linestyles},
                    'linewidth':{'type':'doublespinbox','default':1.0,'range':(0,20),'interval':.2,'label':'line width'},
                    'ms':{'type':'spinbox','default':30,'range':(1,120),'interval':1,'label':'marker size'},
                    'colormap':{'type':'combobox','default':'CMRmap','items':colormaps},
                    'alpha':{'type':'doublespinbox','default':0.9,'range':(.1,1),'interval':.1,'label':'alpha'},
                    'style':{'type':'combobox','default':core.PLOTSTYLE,'items': style_list},
                    'title':{'type':'entry','default':''},
                    'xlabel':{'type':'entry','default':''},
                    'ylabel':{'type':'entry','default':''},
                    'showxlabels':{'type':'checkbox','default':1,'label':'x tick labels'},
                    'showylabels':{'type':'checkbox','default':1,'label':'y tick labels'},
                    'font':{'type':'font','default':defaultfont,'items':fonts},
                    'fontsize':{'type':'spinbox','default':10,'range':(4,50),'label':'font size'},
                    }
        return

class PlotViewer(QWidget):
    """matplotlib plots widget"""
    def __init__(self, parent=None):

        super(PlotViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 900, 600))
        self.main = QSplitter()
        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.main)
        self.create_figure()
        self.create_controls()
        self.setWindowTitle('plots')
        return

    def create_figure(self, fig=None):
        """Create canvas and figure"""

        import matplotlib.pyplot as plt
        #ax.plot(range(10))
        if fig == None:
            fig, ax = plt.subplots(1,1, figsize=(7,5), dpi=120)
            self.ax = ax
        if hasattr(self, 'canvas'):
            self.layout().removeWidget(self.canvas)
        canvas = FigureCanvas(fig)
        left = QWidget()
        self.main.addWidget(left)
        l = QVBoxLayout()
        left.setLayout(l)
        l.addWidget(canvas)
        self.toolbar = NavigationToolbar(canvas, self)
        l.addWidget(self.toolbar)
        self.fig = fig
        self.canvas = canvas
        iconfile = os.path.join(iconpath,'reduce')
        a = QAction(QIcon(iconfile), "Reduce elements",  self)
        a.triggered.connect(lambda: self.zoom(zoomin=False))
        self.toolbar.addAction(a)
        iconfile = os.path.join(iconpath,'enlarge')
        a = QAction(QIcon(iconfile), "Enlarge elements",  self)
        a.triggered.connect(lambda: self.zoom(zoomin=True))
        self.toolbar.addAction(a)
        return

    def create_controls(self):
        """Make widgets for options"""

        self.opts = PlotOptions()
        right = QWidget()
        self.main.addWidget(right)
        w = self.opts.showDialog(self, style=None, section_wrap=1)
        right.setMaximumWidth(220)
        l = QVBoxLayout()
        right.setLayout(l)
        l.addWidget(w)
        btn = QPushButton('Refresh')
        btn.clicked.connect(self.replot)
        l.addWidget(btn)
        btn = QPushButton('Scratchpad')
        iconfile = os.path.join(iconpath,'scratchpad')
        btn.setIcon(QIcon(iconfile))
        btn.clicked.connect(self.save_to_scratchpad)
        l.addWidget(btn)
        self.main.addWidget(right)
        return

    def plot(self, data, kind=None):
        """Do plot"""

        self.opts.applyOptions()
        kwds = self.opts.kwds
        if kind == None and kwds['kind'] != '':
            #overrides kind argument
            kind = kwds['kind']
        else:
            #set widget to current plot kind
            self.opts.setWidgetValue('kind', kind)
            pass
        title = kwds['title']
        #xlabel = kwds['xlabel']
        font = kwds['font']
        fontsize = kwds['fontsize']
        alpha = kwds['alpha']
        cmap = kwds['colormap']
        grid = kwds['grid']
        marker = kwds['marker']
        ms = kwds['ms']
        ls = kwds['linestyle']
        lw = kwds['linewidth']
        subplots = kwds['subplots']
        sharex = kwds['sharex']
        self.style = kwds['style']
        self.set_style()
        self.data = data

        d = data._get_numeric_data()
        xcol = d.columns[0]
        ycols = d.columns[1:]

        nrows = int(round(np.sqrt(len(data.columns)),0))
        layout = (nrows,-1)
        self.clear()
        fig = self.fig
        ax = self.ax
        plt.rc("font", family=kwds['font'], size=fontsize)
        if kind == 'line':
            d.plot(kind='line',ax=ax, cmap=cmap, grid=grid, alpha=alpha, linewidth=lw,
                    linestyle=ls, fontsize=fontsize, subplots=subplots)
        elif kind == 'bar':
            d.plot(kind='bar',ax=ax, cmap=cmap, grid=grid, alpha=alpha, linewidth=lw,
                    fontsize=fontsize, subplots=subplots, sharex=sharex)
        elif kind == 'barh':
            d.plot(kind='barh',ax=ax, cmap=cmap, grid=grid, alpha=alpha, linewidth=lw,
                    fontsize=fontsize, subplots=subplots)
        elif kind == 'hist':
            d.plot(kind='hist',subplots=True,ax=ax,bins=kwds['bins'], linewidth=lw,
                    cmap=cmap, grid=grid, alpha=alpha, fontsize=fontsize)
        elif kind == 'scatter':
            d=d.dropna()
            if len(d.columns)==1:
                x = d.index
                y = d[xcol]
            else:
                x = d[xcol]
                y = d[ycols[0]]
            ax.scatter(x, y, marker=marker, alpha=alpha, linewidth=lw, s=ms) # facecolor=clr,)

        elif kind == 'heatmap':
            #ax.imshow(d, cmap=cmap)
            self.heatmap(d, ax, cmap=cmap, alpha=alpha)
        elif kind == 'pie':
            d.plot(kind='pie',subplots=True,legend=False,layout=layout,ax=ax)
        elif kind == 'box':
            d.boxplot(ax=ax, grid=grid)

        fig.suptitle(title, font=font)
        self.set_axis_labels(ax, kwds)
        plt.tight_layout()
        self.redraw()
        return

    def replot(self):
        """Update current plot"""

        self.plot(self.data, kind=None)
        return

    def refresh(self):
        """Update current plot"""

        self.plot(self.data, kind=None)
        return

    def set_axis_labels(self, ax, kwds):
        """Set a plots axis labels"""

        if kwds['xlabel'] != '':
            ax.set_xlabel(kwds['xlabel'])
        if kwds['ylabel'] != '':
            ax.set_ylabel(kwds['ylabel'])
        if kwds['showxlabels'] == False:
            ax.set_xticklabels([])
        if kwds['showylabels'] == False:
            ax.set_yticklabels([])
        #if kwds['rotx'] != 0:
        #    for tick in ax.get_xticklabels():
        #        tick.set_rotation(kwds['rotx'])
        return

    def heatmap(self, df, ax, cmap='Blues', alpha=0.9, lw=1,
                colorbar=True, cscale=None, grid=False):
        """Plot heatmap"""

        X = df._get_numeric_data()
        clr='black'
        if lw==0:
            clr=None
            lw=None
        if cscale == 'log':
            norm=mpl.colors.LogNorm()
        else:
            norm=None
        hm = ax.pcolor(X, cmap=cmap,linewidth=lw,alpha=alpha,norm=norm)
        if colorbar == True:
            self.fig.colorbar(hm, ax=ax)
        ax.set_xticks(np.arange(0.5, len(X.columns)))
        ax.set_yticks(np.arange(0.5, len(X.index)))
        ax.set_xticklabels(X.columns, minor=False)
        ax.set_yticklabels(X.index, minor=False)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylim(0, len(X.index))
        ax.grid(grid)
        return

    def violinplot(self, df, ax, kwds):
        """violin plot"""

        data=[]
        clrs=[]
        df = df._get_numeric_data()
        cols = len(df.columns)
        cmap = plt.cm.get_cmap(kwds['colormap'])
        for i,d in enumerate(df):
            clrs.append(cmap(float(i)/cols))
            data.append(df[d].values)
        lw = kwds['linewidth']
        alpha = kwds['alpha']
        parts = ax.violinplot(data, showextrema=False, showmeans=True)
        i=0
        for pc in parts['bodies']:
            pc.set_facecolor(clrs[i])
            pc.set_edgecolor('black')
            pc.set_alpha(alpha)
            pc.set_linewidth(lw)
            i+=1
        labels = df.columns
        ax.set_xticks(np.arange(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        return

    def scatter(self, df, ax, axes_layout='single', alpha=0.8, marker='o', color=None, **kwds):
        """A custom scatter plot rather than the pandas one. By default this
        plots the first column selected versus the others"""

        if len(df.columns)<2:
            return
        data = df
        df = df.copy()._get_numeric_data()
        cols = list(df.columns)
        x = df[cols[0]]
        s=1
        cmap = plt.cm.get_cmap(kwds['colormap'])
        lw = kwds['linewidth']
        grid = kwds['grid']
        bw = kwds['bw']

        if cscale == 'log':
            norm = mpl.colors.LogNorm()
        else:
            norm = None

        plots = len(cols)
        if marker == '':
            marker = 'o'
        if axes_layout == 'multiple':
            size = plots-1
            nrows = int(round(np.sqrt(size),0))
            ncols = int(np.ceil(size/nrows))
            self.fig.clear()
        if c is not None:
            colormap = kwds['colormap']
        else:
            colormap = None
            c=None

        #print (kwds)
        labelcol = kwds['labelcol']
        pointsizes = kwds['pointsizes']
        handles = []
        for i in range(s,plots):
            y = df[cols[i]]
            ec = 'black'
            if bw == True:
                clr = 'white'
                colormap = None
            else:
                clr = cmap(float(i)/(plots))
            if colormap != None:
                clr=None
            if marker in ['x','+'] and bw == False:
                ec = clr

            if axes_layout == 'multiple':
                ax = self.fig.add_subplot(nrows,ncols,i)
            if pointsizes != '' and pointsizes in df.columns:
                ms = df[pointsizes]
                s=kwds['ms']
                getsizes = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*s)**2.3
                ms = getsizes(ms)
                #print (ms)
            else:
                ms = kwds['ms'] * 12
            sc = ax.scatter(x, y, marker=marker, alpha=alpha, linewidth=lw, c=c,
                       s=ms, edgecolors=ec, facecolor=clr, cmap=colormap,
                       norm=norm, label=cols[i], picker=True)

            if kwds['logx'] == 1:
                ax.set_xscale('log')
            if kwds['logy'] == 1:
                ax.set_yscale('log')

            #create proxy artist for markers so we can return these handles if needed
            mkr = Line2D([0], [0], marker=marker, alpha=alpha, ms=10, markerfacecolor=c,
                        markeredgewidth=lw, markeredgecolor=ec, linewidth=0)
            handles.append(mkr)
            ax.set_xlabel(cols[0])
            if grid == 1:
                ax.grid(True, linestyle='--')
            else:
                ax.grid(False)
            if axes_layout == 'multiple':
                ax.set_title(cols[i])
            if colormap is not None and kwds['colorbar'] == True:
                self.fig.colorbar(scplt, ax=ax)

            if labelcol != '':
                if not labelcol in data.columns:
                    self.showWarning('label column %s not in selected data' %labelcol)
                elif len(data)<1500:
                    for i, r in data.iterrows():
                        txt = r[labelcol]
                        if pd.isnull(txt) is True:
                            continue
                        ax.annotate(txt, (x[i],y[i]), xycoords='data',
                                    xytext=(5, 5), textcoords='offset points',)

        if kwds['legend'] == 1 and axes_layout == 'single':
            leg = ax.legend(cols[1:])
            #leg.set_draggable(state=True)

        return ax, handles

    def set_figure(self, fig):
        """Set the figure if we have plotted elsewhere"""

        self.clear()
        self.create_figure(fig)
        #self.ax = fig.ax
        self.canvas.draw()
        return

    def clear(self):
        """Clear plot"""

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.canvas.draw()
        return

    def set_style(self):
        """Apply style"""

        if self.style == None or self.style == '':
            mpl.rcParams.update(mpl.rcParamsDefault)
        else:
            plt.style.use(self.style)
        return

    def redraw(self):
        self.canvas.draw()

    def zoom(self, zoomin=True):
        """Zoom in/out to plot by changing size of elements"""

        if zoomin == False:
            val=-1.0
        else:
            val=1.0
        if len(self.opts.kwds) == 0:
            return
        self.opts.increment('linewidth',val/5)
        self.opts.increment('ms',val*2)
        self.opts.increment('fontsize',val)
        self.replot()
        return

    def save_to_scratchpad(self, label=None):
        """Send fig to app scratchpad"""

        print (self.app)
        if not hasattr(self, 'app') or self.app == None:
            return
        name = 'plot'
        if label == None or label is False:
            t = time.strftime("%H:%M:%S")
            label = name+'-'+t
        #get the current figure and make a copy of it by using pickle
        p = pickle.dumps(self.fig)
        fig = pickle.loads(p)
        items = self.app.scratch_items
        items[label] = fig
        if hasattr(self.app, 'scratchpad'):
            self.app.scratchpad.update(items)
        return

class BrowserViewer(QDialog):
    """Browser widget"""
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

        self.browser = QWebEngineView()
        self.browser.urlChanged.connect(self.update_urlbar)
        vbox = QVBoxLayout()
        self.setLayout(vbox)

        # create QToolBar for navigation
        navtb = QToolBar("Navigation")
        vbox.addWidget(navtb)
        # add actions to the tool bar
        iconfile = os.path.join(iconpath,'arrow-left.png')
        icon = QIcon(iconfile)
        back_btn = QAction(icon, "Back", self)
        back_btn.setStatusTip("Back to previous page")
        back_btn.triggered.connect(self.browser.back)
        navtb.addAction(back_btn)
        iconfile = os.path.join(iconpath,'arrow-right.png')
        icon = QIcon(iconfile)
        next_btn = QAction(icon, "Forward", self)
        next_btn.setStatusTip("Forward to next page")
        navtb.addAction(next_btn)
        next_btn.triggered.connect(self.browser.forward)
        iconfile = os.path.join(iconpath,'reload.png')
        icon = QIcon(iconfile)
        reload_btn = QAction(icon, "Reload", self)
        reload_btn.setStatusTip("Reload page")
        reload_btn.triggered.connect(self.browser.reload)
        navtb.addAction(reload_btn)
        navtb.addSeparator()
        # creating a line edit for the url
        self.urlbar = QLineEdit()
        # adding action when return key is pressed
        self.urlbar.returnPressed.connect(self.navigate_to_url)
        navtb.addWidget(self.urlbar)

        vbox.addWidget(self.browser)
        self.browser.setMinimumHeight(500)

        toolswidget = QWidget()
        toolswidget.setMaximumHeight(40)
        vbox.addWidget(toolswidget)
        l = QVBoxLayout(toolswidget)

        self.zoomslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(3)
        w.setMinimum(2)
        w.setMaximum(20)
        w.setValue(10)
        l.addWidget(w)
        w.valueChanged.connect(self.zoom)
        return

    def load_page(self, url):
        self.browser.setUrl(url)

    def navigate_to_url(self):
        """method called by the line edit when return key is pressed
          getting url and converting it to QUrl object"""
        q = QUrl(self.urlbar.text())
        if q.scheme() == "":
            q.setScheme("http")
        # set the url to the browser
        self.browser.setUrl(q)

    def update_urlbar(self, q):
        """method for updating url
           this method is called by the QWebEngineView object
        """
        # setting text to the url bar
        self.urlbar.setText(q.toString())

        # setting cursor position of the url bar
        self.urlbar.setCursorPosition(0)

    def zoom(self):
        zoom = self.zoomslider.value()/10
        self.browser.setZoomFactor(zoom)

class TreeViewer(QWidget):
    def __init__(self, parent=None):

        super(TreeViewer, self).__init__()
        self.parent = parent
        self.setMinimumSize(400,300)
        self.setGeometry(QtCore.QRect(200, 200, 600, 600))
        self.setWindowTitle("Tree View")
        self.add_widgets()
        self.width=600
        self.height=500
        self.treefile = None
        return

    def add_widgets(self):

        layout = QVBoxLayout(self)
        self.browser = QWebEngineView()
        layout.addWidget(self.browser)
        toolswidget = QWidget()
        layout.addWidget(toolswidget)
        toolswidget.setMaximumHeight(120)
        l2 = QVBoxLayout(toolswidget)

        self.zoomslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(3)
        w.setMinimum(2)
        w.setMaximum(20)
        w.setValue(10)
        l2.addWidget(w)
        w.valueChanged.connect(self.zoom)

        w = QLabel('vscale')
        l2.addWidget(w)
        self.heightslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(10)
        w.setMinimum(20)
        w.setMaximum(1200)
        w.setValue(500)
        l2.addWidget(w)
        w.valueChanged.connect(self.vscale)

        w = QLabel('hscale')
        l2.addWidget(w)
        self.widthslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(10)
        w.setMinimum(20)
        w.setMaximum(2000)
        w.setValue(600)
        l2.addWidget(w)
        w.valueChanged.connect(self.hscale)
        return

    def draw(self, treefile, df, **kwargs):
        """Show phylogeny"""

        import toyplot
        from . import trees
        self.treefile = treefile
        self.df = df
        canvas = trees.draw_tree(treefile,  df,
                             width=self.width, height=self.height, **kwargs)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            self.browser.setHtml(html)
        return

    def update(self):

        if self.treefile == None:
            return
        self.draw(self.treefile, self.df, tiplabelcol='species')
        return

    def zoom(self):
        zoom = self.zoomslider.value()/10
        self.browser.setZoomFactor(zoom)
        return

    def hscale(self):
        self.width = self.widthslider.value()
        self.update()
        return

    def vscale(self):
        self.height = self.heightslider.value()
        self.update()
        return

class TableViewer(QDialog):
    """View row of data in table"""
    def __init__(self, parent=None, dataframe=None, **kwargs):
        super(TableViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 600, 600))
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.table = tables.DataFrameTable(self, dataframe=dataframe, **kwargs)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.grid.addWidget(self.table)
        return

class AlignmentWidget(QWidget):
    """Widget for showing sequence alignments"""
    def __init__(self, parent=None):
        super(AlignmentWidget, self).__init__(parent)
        l = QHBoxLayout(self)
        self.setLayout(l)
        self.m = QSplitter(self)
        l.addWidget(self.m)
        self.left = PlainTextEditor(self.m, readOnly=True)
        self.right = PlainTextEditor(self.m, readOnly=True)
        self.left.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.right.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.m.setSizes([200,300])
        self.m.setStretchFactor(1,2)
        return

class SequencesViewer(QWidget):
    """Viewer for sequences and alignments"""

    def __init__(self, parent=None, filename=None):
        super(SequencesViewer, self).__init__(parent)

        #self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setGeometry(QtCore.QRect(200, 200, 1000, 600))
        self.setMinimumHeight(150)
        self.recs = None
        self.aln = None
        self.add_widgets()
        self.show()
        return

    def add_widgets(self):
        """Add widgets"""

        #self.main = QWidget(self)
        #self.setCentralWidget(self.main)
        l = QHBoxLayout(self)
        self.setLayout(l)
        self.tabs = QTabWidget(self)
        l.addWidget(self.tabs)

        self.ed = ed = PlainTextEditor(self, readOnly=True)
        self.ed.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.tabs.addTab(self.ed, 'fasta')

        self.alnview = AlignmentWidget(self)
        self.tabs.addTab(self.alnview, 'alignment')

        sidebar = QWidget()
        sidebar.setFixedWidth(180)
        l.addWidget(sidebar)
        l2 = QVBoxLayout(sidebar)
        l2.setSpacing(5)
        l2.setAlignment(QtCore.Qt.AlignTop)
        btn = create_button(self, None, self.zoom_out,
                                    'zoom-out', core.ICONSIZE, 'zoom out')
        l2.addWidget(btn)
        btn = create_button(self, None, self.zoom_in,
                                    'zoom-in', core.ICONSIZE, 'zoom in')
        l2.addWidget(btn)
        lbl = QLabel('Format')
        l2.addWidget(lbl)
        w = QComboBox()
        w.addItems(['no color','color by residue','color by difference'])
        w.setCurrentIndex(1)
        w.activated.connect(self.show_alignment)
        self.formatchoice = w
        l2.addWidget(w)
        lbl = QLabel('Set Reference')
        l2.addWidget(lbl)
        w = QComboBox()
        w.activated.connect(self.set_reference)
        self.referencechoice = w
        l2.addWidget(w)
        lbl = QLabel('Aligner')
        l2.addWidget(lbl)
        w = QComboBox()
        w.setCurrentIndex(1)
        w.addItems(['clustal','muscle'])
        self.alignerchoice = w
        l2.addWidget(w)
        #self.create_menu()
        return

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('&File', self)
        self.menuBar().addMenu(self.file_menu)
        self.file_menu.addAction('&Load Fasta File', self.load_fasta,
                QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.file_menu.addAction('&Save Alignment', self.save_alignment,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        return

    def scroll_top(self):
        vScrollBar = self.ed.verticalScrollBar()
        vScrollBar.triggerAction(QScrollBar.SliderToMinimum)
        return

    def zoom_out(self):
        self.ed.zoom(-1)
        self.alnview.left.zoom(-1)
        self.alnview.right.zoom(-1)
        return

    def zoom_in(self):
        self.ed.zoom(1)
        self.alnview.left.zoom(1)
        self.alnview.right.zoom(1)
        return

    def load_fasta(self, filename=None):
        """Load fasta file"""

        if filename == None:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                            filter="Fasta Files(*.fa *.fna *.fasta);;All Files(*.*)")
        if not filename:
            return
        from Bio import AlignIO, SeqIO
        recs = list(SeqIO.parse(filename, 'fasta'))
        self.load_records(recs)
        return

    def load_records(self, recs):
        """Load seqrecords list"""

        from Bio import AlignIO, SeqIO
        self.recs = recs
        self.reference = self.recs[0]
        rdict = SeqIO.to_dict(recs)
        self.show_fasta()
        self.show_alignment()
        self.referencechoice.addItems(list(rdict.keys()))
        return

    def load_alignment(self, aln_file):
        """Load alignment directly from a fasta file"""

        from Bio import AlignIO, SeqIO
        self.aln = AlignIO.read(aln_file, 'fasta')
        recs = tools.seqrecords_from_alignment(self.aln)
        self.load_records(recs)
        #self.show_alignment()
        return

    def set_reference(self):
        ref = self.referencechoice.currentText()
        return

    def show_fasta(self):
        """Show records as fasta"""

        recs = self.recs
        if recs == None:
            return
        self.ed.clear()
        for rec in recs:
            s = rec.format('fasta')
            self.ed.insertPlainText(s)
        self.scroll_top()
        return

    def align(self):
        """Align current sequences"""

        from Bio import AlignIO, SeqIO
        if self.aln == None:
            outfile = 'temp.fa'
            SeqIO.write(self.recs, outfile, 'fasta')
            self.aln = tools.clustal_alignment(outfile)
        return

    def show_alignment(self):

        format = self.formatchoice.currentText()
        self.draw_alignment(format)
        return

    def draw_alignment(self, format='color by residue'):
        """Show alignment with colored columns"""

        left = self.alnview.left
        right = self.alnview.right
        chunks=0
        offset=0
        diff=False
        self.align()
        aln = self.aln
        left.clear()
        right.clear()
        self.scroll_top()

        #colors = tools.get_protein_colors()
        colors = {
            'A': 'lightgreen',
            'T': 'lightblue',
            'C': 'lightyellow',
            'G': 'lightcoral',
            '-': 'white',
            'N': 'gray'
        }

        format = QTextCharFormat()
        format.setBackground(QBrush(QColor('white')))
        cursor = right.textCursor()

        ref = aln[0]
        l = len(aln[0])
        n=60
        s=[]
        if chunks > 0:
            chunks = [(i,i+n) for i in range(0, l, n)]
        else:
            chunks = [(0,l)]
        for c in chunks:
            start,end = c
            lbls = np.arange(start+1,end+1,10)-offset
            head = ''.join([('%-10s' %i) for i in lbls])
            cursor.insertText(head)
            right.insertPlainText('\n')
            left.appendPlainText(' ')
            for a in aln:
                name = a.id
                seq = a.seq[start:end].upper()
                left.appendPlainText(name)
                line = ''
                for i in seq:
                    if i in colors:
                        c = colors[i]
                    else:
                        c = 'white'
                    line += '<span style="background-color:%s;">%s</span>' %(c,i)
                cursor.insertHtml(line)
                right.insertPlainText('\n')
        return

    def save_alignment(self):

        filters = "clustal files (*.aln);;All files (*.*)"
        filename, _ = QFileDialog.getSaveFileName(self,"Save Alignment",
                                                  "",filters)
        if not filename:
            return
        SeqIO.write(self.aln,filename,format='clustal')
        return

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

        self.startw = QLineEdit()
        bl.addWidget(self.startw)
        self.startw.setMaximumWidth(100)
        self.startw.returnPressed.connect(self.goto)

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
        self.length = length = plotting.get_fasta_length(ref_file)
        if self.gb_file != None:
            df = tools.genbank_to_dataframe(gb_file)
            if 'locus_tag' in df.columns:
                df.loc[df["gene"].isnull(),'gene'] = df.locus_tag
            if 'gene' in df.columns:
                genes = df.gene.dropna().unique()
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
        sl.setTickInterval(int(length/20))
        return

    def set_chrom(self, chrom):
        """Set the selected record which also updates the plot"""

        index = self.chromselect.findText(chrom)
        self.chromselect.setCurrentIndex(index)
        return

    def update_chrom(self, chrom=None):
        """Update after chromosome selection changed"""

        #self.chrom = chrom
        recname = self.chromselect.currentText()
        length = self.length = plotting.get_fasta_length(self.ref_file, key=chrom)
        sl = self.slider
        sl.setMinimum(1)
        sl.setMaximum(length)
        sl.setTickInterval(int(length/20))
        self.redraw()
        return

    def goto(self):

        loc = int(self.startw.text())
        self.redraw(loc)
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

class ScratchPad(QWidget):
    """Temporary storage widget for plots and other items.
    Currently supports storing text, mpl figures and dataframes"""
    def __init__(self, parent=None):
        super(ScratchPad, self).__init__(parent)
        self.parent = parent
        self.setMinimumSize(400,300)
        self.setGeometry(QtCore.QRect(300, 200, 800, 600))
        self.setWindowTitle("Scratchpad")
        self.createWidgets()
        sizepolicy = QSizePolicy()
        self.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))
        #dict to store objects, these should be serialisable
        self.items = {}
        return

    def createWidgets(self):
        """Create widgets. Plot on left and dock for tools on right."""

        self.main = QTabWidget(self)
        self.main.setTabsClosable(True)
        self.main.tabCloseRequested.connect(lambda index: self.remove(index))
        layout = QVBoxLayout(self)
        toolbar = QToolBar("toolbar")
        layout.addWidget(toolbar)
        items = { 'new text':{'action':self.newText,'file':'document-new'},
                  'save': {'action':self.save,'file':'save'},
                  'save all': {'action':self.saveAll,'file':'save-all'},
                  'clear': {'action':self.clear,'file':'clear'}
                    }
        for i in items:
            if 'file' in items[i]:
                iconfile = os.path.join(iconpath,items[i]['file']+'.png')
                icon = QIcon(iconfile)
            else:
                icon = QIcon.fromTheme(items[i]['icon'])
            btn = QAction(icon, i, self)
            btn.triggered.connect(items[i]['action'])
            toolbar.addAction(btn)
        layout.addWidget(self.main)
        return

    def update(self, items):
        """Display a dict of stored objects"""

        self.main.clear()
        for name in items:
            obj = items[name]
            #print (name,type(obj))
            if type(obj) is str:
                te = dialogs.PlainTextEditor()
                te.setPlainText(obj)
                self.main.addTab(te, name)
            elif type(obj) is pd.DataFrame:
                tw = core.DataFrameTable(self.main, dataframe=obj)
                self.main.addTab(tw, name)
            else:
                pw = PlotWidget(self.main)
                self.main.addTab(pw, name)
                pw.figure = obj
                pw.draw()
                plt.tight_layout()
        self.items = items
        return

    def remove(self, idx):
        """Remove selected tab and item widget"""

        index = self.main.currentIndex()
        name = self.main.tabText(index)
        del self.items[name]
        self.main.removeTab(index)
        return

    def save(self):
        """Save selected item"""

        index = self.main.currentIndex()
        name = self.main.tabText(index)
        suff = "PNG files (*.png);;JPG files (*.jpg);;PDF files (*.pdf);;All files (*.*)"
        filename, _ = QFileDialog.getSaveFileName(self, "Save Figure", name, suff)
        if not filename:
            return

        fig = self.items[name]
        fig.savefig(filename+'.png', dpi=core.DPI)
        return

    def saveAll(self):
        """Save all figures in a folder"""

        dir =  QFileDialog.getExistingDirectory(self, "Save Folder",
                                             homepath, QFileDialog.ShowDirsOnly)
        if not dir:
            return
        for name in self.items:
            fig = self.items[name]
            fig.savefig(os.path.join(dir,name+'.png'), dpi=core.DPI)
        return

    def clear(self):
        """Clear plots"""

        self.items.clear()
        self.main.clear()
        return

    def newText(self):
        """Add a text editor"""

        name, ok = QInputDialog.getText(self, 'Name', 'Name:',
                    QLineEdit.Normal, '')
        if ok:
            tw = dialogs.PlainTextEditor()
            self.main.addTab(tw, name)
            self.items[name] = tw.toPlainText()
        return

    def closeEvent(self, event):
        """Close"""

        for idx in range(self.main.count()):
            name = self.main.tabText(idx)
            #print (name)
            w = self.main.widget(idx)
            #print (w)
            if type(w) == PlainTextEditor:
                self.items[name] = w.toPlainText()
        return
