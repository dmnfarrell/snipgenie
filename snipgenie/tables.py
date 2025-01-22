#!/usr/bin/env python

"""
    dataframe table widget and sub-classes.
    Created Nov 2019
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

import sys,os,platform
import pandas as pd
import numpy as np
from .qt import *
from pandas.api.types import is_datetime64_any_dtype as is_datetime
import pylab as plt
import matplotlib as mpl
from . import widgets, core, tools, plotting

class ColumnHeader(QHeaderView):
    def __init__(self):
        super(QHeaderView, self).__init__()
        return

class DataFrameWidget(QWidget):
    """Widget containing a tableview and statusbar"""
    def __init__(self, parent=None, table=None, statusbar=True, toolbar=False, app=None, **kwargs):
        """
        Widget containing a dataframetable - allows us to pass any table subclass
        """

        super(DataFrameWidget, self).__init__()
        l = self.layout = QGridLayout()
        l.setSpacing(2)
        self.setLayout(self.layout)
        #self.plotview = widgets.PlotViewer()
        if table == None:
            self.table = DataFrameTable(self, dataframe=pd.DataFrame(), app=app)
        else:
            self.table = table
            table.parent = self
        l.addWidget(self.table, 1, 1)
        if toolbar == True:
            self.createToolbar()
        if statusbar == True:
            self.statusBar()

        self.table.model.dataChanged.connect(self.stateChanged)
        return

    #@Slot('QModelIndex','QModelIndex','int')
    def stateChanged(self, idx, idx2):
        """Run whenever table model is changed"""

        if hasattr(self, 'pf') and self.pf is not None:
            self.pf.updateData()

    def statusBar(self):
        """Status bar at bottom"""

        w = self.statusbar = QWidget(self)
        l = QHBoxLayout(w)
        w.setMaximumHeight(30)
        self.size_label = QLabel("")
        l.addWidget(self.size_label, 1)
        w.setStyleSheet('color: #1a216c; font-size:12px')
        self.layout.addWidget(w, 2, 1)
        self.updateStatusBar()
        return

    def updateStatusBar(self):
        """Update the table details in the status bar"""

        if not hasattr(self, 'size_label'):
            return
        df = self.table.model.df
        #meminfo = self.table.getMemory()
        s = '{r} samples x {c} columns'.format(r=len(df), c=len(df.columns))
        self.size_label.setText(s)
        return

    def createToolbar(self):

        self.setLayout(self.layout)
        items = {
                 'copy': {'action': self.copy,'file':'copy','shortcut':'Ctrl+C'},
                 'line': {'action': lambda: self.table.plot(kind='line'),'file':'plot_line'},
                 'bar': {'action': lambda: self.table.plot(kind='bar'),'file':'plot_bar'},
                 'barh': {'action': lambda: self.table.plot(kind='barh'),'file':'plot_barh'},
                 'hist': {'action': lambda: self.table.plot(kind='hist'),'file':'plot_hist'},
                 'scatter': {'action': lambda: self.table.plot(kind='scatter'),'file':'plot_scatter'},
                 'heatmap': {'action': lambda: self.table.plot(kind='heatmap'), 'file':'plot_heatmap'},
                 'pie': {'action': lambda: self.table.plot(kind='pie'),'file':'plot_pie'},
                 #'filter':{'action':sel-f.filter,'file':'table-filter'}
                 }

        self.toolbar = toolbar = QToolBar("Toolbar")
        toolbar.setIconSize(QtCore.QSize(core.ICONSIZE-4, core.ICONSIZE-4))
        toolbar.setOrientation(QtCore.Qt.Vertical)
        widgets.addToolBarItems(toolbar, self, items)
        self.layout.addWidget(toolbar,1,2)
        return

    def refresh(self):

        self.table.refresh()
        return

    def copy(self):
        """Copy to clipboard"""

        #check size of dataframe
        m = self.table.model.df.memory_usage(deep=True).sum()
        if m>1e8:
            answer = QMessageBox.question(self, 'Copy?',
                             'This data may be too large to copy. Are you sure?', QMessageBox.Yes, QMessageBox.No)
            if not answer:
                return
        df = self.table.getSelectedDataFrame()
        df.to_clipboard()
        return

    def filter(self):
        """Show filter dialog"""

        return

class HeaderView(QHeaderView):
    """"
    Column header class.
    """
    def __init__(self, parent):
        super(HeaderView, self).__init__(QtCore.Qt.Horizontal, parent)
        '''self.setStyleSheet(
            "QHeaderView::section{background-color: #ffffff; "
            "font-weight: bold; "
            "border-bottom: 1px solid gray;}")'''

        self.setDefaultAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.Alignment(QtCore.Qt.TextWordWrap))
        sizePol = QSizePolicy()
        sizePol.setVerticalPolicy(QSizePolicy.Maximum)
        sizePol.setHorizontalPolicy(QSizePolicy.Maximum)
        self.setSizePolicy(sizePol)
        self.MAX_HEIGHT = 240
        self.setMinimumHeight(26)
        self.setMaximumHeight(self.MAX_HEIGHT )
        self.setSectionsClickable(True)
        self.setSelectionBehavior(QTableView.SelectColumns)
        self.setStretchLastSection(False)
        return

    def sectionSizeFromContents(self, logicalIndex):
        """Get section size from contents"""

        text = self.model().headerData(logicalIndex, self.orientation(), QtCore.Qt.DisplayRole)
        alignment = self.defaultAlignment()
        metrics = QFontMetrics(self.fontMetrics())
        width = metrics.boundingRect(QtCore.QRect(), alignment, text).width()

        heights = []
        for i in range(self.count()):
            text = self.model().headerData(i, self.orientation(),QtCore.Qt.DisplayRole)
            size = self.sectionSize(i)
            rect = QtCore.QRect(0, 0, size, self.MAX_HEIGHT)
            heights.append(metrics.boundingRect(rect, alignment, text).height())
        height = sorted(heights)[-1] + 5
        return QtCore.QSize(width, height)

class DataFrameTable(QTableView):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, plotter=None,
                    font=core.FONT, fontsize=10):

        QTableView.__init__(self)
        self.parent = parent
        self.clicked.connect(self.showSelection)

        vh = self.verticalHeader()
        vh.setVisible(True)
        vh.setDefaultSectionSize(28)
        vh.setMinimumWidth(20)
        vh.setMaximumWidth(500)

        self.headerview = HeaderView(self)
        self.setHorizontalHeader(self.headerview)
        hh = self.horizontalHeader()
        hh.setVisible(True)
        hh.setSectionsMovable(True)
        hh.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        hh.customContextMenuRequested.connect(self.columnHeaderMenu)
        hh.setSelectionBehavior(QTableView.SelectColumns)
        hh.setSelectionMode(QAbstractItemView.ExtendedSelection)
        #hh.sectionClicked.connect(self.columnClicked)
        hh.setSectionsClickable(True)

        self.setDragEnabled(True)
        self.viewport().setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.resizeColumnsToContents()
        self.setCornerButtonEnabled(True)

        self.fontname = font
        self.fontsize = fontsize
        self.updateFont()
        tm = DataFrameModel(dataframe)
        self.setModel(tm)
        self.model = tm
        self.setWordWrap(False)
        self.setCornerButtonEnabled(True)
        if plotter == None:
            self.plotview = widgets.PlotViewer()
        else:
            self.plotview = plotter
        self.plotview.app = app
        return

    def getMaxHeight(self):
        h = max([len(str(c))*10 for c in self.model.df.columns])
        #if h>150: h=150
        return h

    def updateFont(self):
        """Update the font"""

        font = self.font = QFont(self.fontname)
        font.setPointSize(int(self.fontsize))
        self.setFont(font)
        self.horizontalHeader().setFont(font)
        self.verticalHeader().setFont(font)
        return

    def setDataFrame(self, df):

        tm = DataFrameModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def getDataFrame(self):
        return self.model.df

    def zoomIn(self, fontsize=None):

        if fontsize == None:
            s = self.font.pointSize()+1
        else:
            s = fontsize
        self.font.setPointSize(s)
        self.setFont(self.font)
        vh = self.verticalHeader()
        h = vh.defaultSectionSize()
        vh.setDefaultSectionSize(h+2)
        hh = self.horizontalHeader()
        w = hh.defaultSectionSize()
        hh.setDefaultSectionSize(w+2)
        return

    def zoomOut(self, fontsize=None):

        if fontsize == None:
            s = self.font.pointSize()-1
        else:
            s = fontsize
        self.font.setPointSize(s)
        self.setFont(self.font)
        vh = self.verticalHeader()
        h = vh.defaultSectionSize()
        vh.setDefaultSectionSize(h-2)
        hh = self.horizontalHeader()
        w = hh.defaultSectionSize()
        hh.setDefaultSectionSize(w-2)
        return

    def importFile(self, filename=None, dialog=False, **kwargs):

        if dialog is True:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Import File",
                                                      "","All Files (*);;Text Files (*.txt);;CSV files (*.csv)",
                                                      options=options)
            df = pd.read_csv(filename, **kwargs)
            self.table.model.df = df
        return

    def info(self):

        buf = io.StringIO()
        self.table.model.df.info(verbose=True,buf=buf,memory_usage=True)
        td = dialogs.TextDialog(self, buf.getvalue(), 'Info')
        return

    def showSelection(self, item):

        cellContent = item.data()
        #print(cellContent)  # test
        row = item.row()
        model = item.model()
        columnsTotal= model.columnCount(None)
        return

    def getSelectedRows(self):

        sm = self.selectionModel()
        rows = [(i.row()) for i in sm.selectedIndexes()]
        rows = list(dict.fromkeys(rows).keys())
        return rows

    def getSelectedIndexes(self):
        """Get selected row indexes"""

        rows = self.getSelectedRows()
        idx = self.model.df.index[rows]
        return idx

    def getSelectedColumns(self):
        """Get selected column indexes"""

        sm = self.selectionModel()
        cols = [(i.column()) for i in sm.selectedIndexes()]
        cols = list(dict.fromkeys(cols).keys())
        return cols

    def getSelectedDataFrame(self):
        """Get selection as a dataframe"""

        df = self.model.df
        sm = self.selectionModel()
        rows = [(i.row()) for i in sm.selectedIndexes()]
        cols = [(i.column()) for i in sm.selectedIndexes()]
        #get unique rows/cols keeping order
        rows = list(dict.fromkeys(rows).keys())
        cols = list(dict.fromkeys(cols).keys())
        return df.iloc[rows,cols]

    def setSelected(self, rows, cols):
        """
        Set selection programmatically from a list of rows and cols.
        https://doc.qt.io/archives/qtjambi-4.5.2_01/com/trolltech/qt/model-view-selection.html
        """

        #print (rows,cols)
        if len(rows)==0 or len(cols)==0:
            return
        topleft = self.model.index(rows[0], cols[0])
        bottomright = self.model.index(rows[-1], cols[-1])
        selection = QtCore.QItemSelection(topleft, bottomright)
        mode = QtCore.QItemSelectionModel.Select
        self.selectionModel().select(selection, mode)
        return

    def getScrollPosition(self):
        """Get current row/col position"""

        hb = self.horizontalScrollBar()
        vb = self.verticalScrollBar()
        return vb.value(),hb.value()

    def setScrollPosition(self, row, col):
        """Move to row/col position"""

        idx = self.model.index(row, col)
        self.scrollTo(idx)
        return

    def handleDoubleClick(self, item):

        cellContent = item.data()
        if item.column() != 0:
            return
        return

    def columnClicked(self, col):

        hheader = self.horizontalHeader()
        df = self.model.df
        return

    def storeCurrent(self):
        """Store current version of the table before a major change is made"""

        self.prevdf = self.model.df.copy()
        return

    def deleteCells(self, rows, cols, answer=None):
        """Clear the cell contents"""

        if answer == None:
            answer = QMessageBox.question(self, 'Delete Cells?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if not answer:
            return
        self.storeCurrent()
        #print (rows, cols)
        self.model.df.iloc[rows,cols] = np.nan
        return

    def editCell(self, item):
        return

    def setRowColor(self, rowIndex, color):
        for j in range(self.columnCount()):
            self.item(rowIndex, j).setBackground(color)

    def columnHeaderMenu(self, pos):

        hheader = self.horizontalHeader()
        idx = hheader.logicalIndexAt(pos)
        column = self.model.df.columns[idx]
        menu = QMenu(self)
        sortAction = menu.addAction("Sort \u2193")
        sortDescAction = menu.addAction("Sort \u2191")
        deleteColumnAction = menu.addAction("Delete Column")
        renameColumnAction = menu.addAction("Rename Column")
        addColumnAction = menu.addAction("Add Column")
        plotAction = menu.addAction("Histogram")
        colorbyAction = menu.addAction("Color By Column")

        action = menu.exec_(self.mapToGlobal(pos))
        if action == sortAction:
            self.sort(idx)
        elif action == sortDescAction:
            self.sort(idx, ascending=False)
        elif action == deleteColumnAction:
            self.deleteColumn(column)
        elif action == renameColumnAction:
            self.renameColumn(column)
        elif action == addColumnAction:
            self.addColumn()
        elif action == plotAction:
            self.plotHist(column)
        elif action == colorbyAction:
            self.colorByColumn(column)
        return

    def keyPressEvent(self, event):

        rows = self.getSelectedRows()
        cols = self.getSelectedColumns()
        if event.key() == QtCore.Qt.Key_Delete:
            self.deleteCells(rows, cols)

    def contextMenuEvent(self, event):
        """Reimplemented to create context menus for cells and empty space."""

        # Determine the logical indices of the cell where click occured
        hheader, vheader = self.horizontalHeader(), self.verticalHeader()
        position = event.globalPos()
        row = vheader.logicalIndexAt(vheader.mapFromGlobal(position))
        column = hheader.logicalIndexAt(hheader.mapFromGlobal(position))
        if row == -1:
            return
        # Show a context menu for empty space at bottom of table...
        self.menu = QMenu(self)
        self.addActions(event, row, column)
        return

    def addActions(self, event, row, column):
        """Actions"""

        menu = self.menu
        copyAction = menu.addAction("Copy")
        exportAction = menu.addAction("Export Table")
        transposeAction = menu.addAction("Transpose")
        action = menu.exec_(event.globalPos())
        if action == copyAction:
            self.copy()
        elif action == exportAction:
            self.exportTable()
        elif action == transposeAction:
            self.transpose()
        return

    def setIndex(self):
        return

    def copy(self):
        """Copy selected cells"""

        df = self.getSelectedDataFrame()
        if len(df.columns==1):
            df.to_clipboard(index=False, header=False)
        else:
            df.to_clipboard(index=False, header=True)
        return

    def refresh(self):
        """Refresh table if dataframe is changed"""

        self.model.beginResetModel()
        index = self.model.index
        try:
            self.model.dataChanged.emit(0,0)
        except:
            self.model.dataChanged.emit(index(0,0),index(0,0))
        self.model.endResetModel()
        if hasattr(self.parent,'statusbar'):
            self.parent.updateStatusBar()
        return

    def importFile(self):
        dialogs.ImportDialog(self)
        return

    def exportTable(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export",
                                                  "","csv files (*.csv);;All files (*.*)")
        if filename:
            self.model.df.to_csv(filename)
        return

    def addColumn(self):
        """Add a  column"""

        df = self.model.df
        name, ok = QInputDialog().getText(self, "Enter Column Name",
                                             "Name:", QLineEdit.Normal)
        if ok and name:
            if name in df.columns:
                return
            df[name] = pd.Series()
            self.refresh()
        return

    def deleteColumn(self, column=None):

        reply = QMessageBox.question(self, 'Delete Columns?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.No:
            return False
        self.model.df = self.model.df.drop(columns=[column])
        self.refresh()
        return

    def deleteRows(self):

        rows = self.getSelectedRows()
        reply = QMessageBox.question(self, 'Delete Rows?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.No:
            return False
        idx = self.model.df.index[rows]
        self.model.df = self.model.df.drop(idx)
        self.refresh()
        return

    def renameColumn(self, column=None):

        name, ok = QInputDialog().getText(self, "Enter New Column Name",
                                             "Name:", QLineEdit.Normal)
        if ok and name:
            self.model.df.rename(columns={column:name},inplace=True)
            self.refresh()
        return

    def sort(self, idx, ascending=True):
        """Sort by selected columns"""

        df = self.model.df
        sel = self.getSelectedColumns()
        if len(sel)>1:
            for i in sel:
                self.model.sort(i, ascending)
        else:
            self.model.sort(idx, ascending)
        return

    def transpose(self):

        self.model.df = self.model.df.T
        self.refresh()
        return

    def plot(self, kind='bar'):
        """Plot table selection"""

        df = self.model.df
        data = self.getSelectedDataFrame()
        self.plotview.plot(data, kind=kind)
        self.plotview.show()
        self.plotview.activateWindow()
        return

    def colorByColumn(self, col):
        """Set colorby column"""

        #cmap = 'Set1'
        df = self.model.df

        colors,colormap = plotting.get_color_mapping(df,col,seed=10, cmap=cmap)
        print (colors)
        self.model.rowcolors = colors
        return

    def getMemory(self):
        """Get memory info as string"""

        m = self.model.df.memory_usage(deep=True).sum()
        if m>1e5:
            m = round(m/1048576,2)
            units='MB'
        else:
            units='Bytes'
        s = "%s %s" %(m,units)
        return s

class DataFrameModel(QtCore.QAbstractTableModel):
    def __init__(self, dataframe=None, *args):
        super(DataFrameModel, self).__init__()
        if dataframe is None:
            self.df = pd.DataFrame()
        else:
            self.df = dataframe
        self.bg = '#F4F4F3'
        self.rowcolors = None
        return

    def update(self, df):
        self.df = df

    def rowCount(self, parent=QtCore.QModelIndex()):

        return len(self.df)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.columns)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        """Edit or display roles. Handles what happens when the Cells
        are edited or what appears in each cell.
        """

        i = index.row()
        j = index.column()
        #print (self.df.dtypes)
        #coltype = self.df.dtypes[j]
        coltype = self.df[self.df.columns[j]].dtype
        isdate = is_datetime(coltype)
        if role == QtCore.Qt.DisplayRole:
            value = self.df.iloc[i, j]
            if isdate:
                return value.strftime(TIMEFORMAT)
            elif type(value) != str:
                if type(value) in [float,np.float64] and np.isnan(value):
                    return ''
                elif type(value) == float:
                    return value
                else:
                    return (str(value))
            else:
                return '{0}'.format(value)
        elif (role == QtCore.Qt.EditRole):
            value = self.df.iloc[i, j]
            #print (coltype)
            #print (value)
            if type(value) is str:
                try:
                    return float(value)
                except:
                    return str(value)

            if np.isnan(value):
                return ''
        elif role == QtCore.Qt.BackgroundRole:
            if self.rowcolors != None:
                clr = self.rowcolors[i]
                #print (clr)
                return QColor(clr)
            else:
                return QColor(self.bg)

    def headerData(self, col, orientation, role=QtCore.Qt.DisplayRole):
        """What's displayed in the headers"""

        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return str(self.df.columns[col])
            if orientation == QtCore.Qt.Vertical:
                value = self.df.index[col]
                if type( self.df.index) == pd.DatetimeIndex:
                    if not value is pd.NaT:
                        try:
                            return value.strftime(TIMEFORMAT)
                        except:
                            return ''
                else:
                    return str(value)
        return None

    def sort(self, idx, ascending=True):
        """Sort table by given column number """

        self.layoutAboutToBeChanged.emit()
        col = self.df.columns[idx]
        self.df = self.df.sort_values(col, ascending=ascending)
        self.layoutChanged.emit()
        return

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        """Set data upon edits"""

        #print (index)
        i = index.row()
        j = index.column()
        curr = self.df.iloc[i,j]
        #print (curr, value)
        self.df.iloc[i,j] = value
        #self.dataChanged.emit()
        return True

    def onDataChanged(self):
        #print (self.df)
        return

    def setColumnColor(self, columnIndex, color):
        for i in range(self.rowCount()):
            self.item(i, columnIndex).setBackground(color)
        return

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled|Qt.ItemIsEditable

class DefaultTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        self.app = app
        self.setWordWrap(False)

class SampleTableModel(DataFrameModel):
    """Samples table model class"""
    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.df = dataframe

    def data(self, index, role):
        """Custom display for sample table"""

        i = index.row()
        j = index.column()
        rowname = self.df.index[i]
        value = self.df.iloc[i, j]
        colname = self.df.columns[j]
        #print (self.df.dtypes)
        coltype = self.df[self.df.columns[j]].dtype
        isdate = is_datetime(coltype)
        color = QColor(self.bg)
        if role == QtCore.Qt.DisplayRole:
            value = self.df.iloc[i, j]
            if isdate:
                return value.strftime(TIMEFORMAT)
            elif type(value) != str:
                if type(value) in [float,np.float64] and np.isnan(value):
                    return ''
                elif type(value) == np.float32:
                    return value
                else:
                    return (str(value))
            else:
                return '{0}'.format(value)
        elif (role == QtCore.Qt.EditRole):
            value = self.df.iloc[i, j]
            #print (coltype)
            #print (value)
            if type(value) is str:
                try:
                    return float(value)
                except:
                    return str(value)

            if np.isnan(value):
                return ''
        elif role == QtCore.Qt.BackgroundRole:
            #color warnings
            if colname == 'meandepth':
                if value < 20:
                    color = QColor('#E57A6D')
                    #color.setAlpha(200)
            elif colname == 'coverage':
                if value < 98:
                    color = QColor('#E08560')
                    #color.setAlpha(200)
            elif colname == 'quality':
                if value == 'pass':
                    color = QColor('#94DA60')
            elif self.rowcolors != None:
                clr = self.rowcolors[i]
                color = QColor(clr)
            elif colname == 'dup':
                if value == True:
                    color = QColor('#f5cf7c')
            else:
                color = QColor(self.bg)
            return color

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class SampleTable(DataFrameTable):
    """
    QTableView for files/samples view.
    """
    def __init__(self, parent=None, app=None, dataframe=None, plotter=None):
        DataFrameTable.__init__(self)
        self.parent = parent
        self.app = app
        self.setWordWrap(False)
        tm = SampleTableModel(dataframe)
        self.setModel(tm)
        self.plotview.app = app
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        if 'sample' in df.columns:
            df = df.set_index('sample', drop=False)
            df.index.name = 'index'
        tm = SampleTableModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def addActions(self, event, row, column):
        """Table actions"""

        menu = self.menu
        detailsAction = menu.addAction("Sample Details")
        fastqqualityAction = menu.addAction("Quality Summary")
        #readlengthsAction = menu.addAction("Read Length Distribution")
        #plotbamAction = menu.addAction("Show Read Alignments")
        mappingstatsAction = menu.addAction("Mapping Statistics")
        #contamaction = menu.addAction('Check Contamination')
        samplereadsaction = menu.addAction('Sample Sequences')
        normalisefastqaction = menu.addAction('Normalise Fastq')
        removeAction = menu.addAction("Remove Selected")
        removebamAction = menu.addAction("Delete Bam Files")
        movefastqAction = menu.addAction("Move Fastq Files")
        copyAction = menu.addAction("Copy")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        # Map the logical row index to a real index for the source model
        #model = self.model
        rows = self.getSelectedRows()
        if action == detailsAction:
            self.app.sample_details(row)
        elif action == fastqqualityAction:
            #print (row)
            self.app.quality_summary(row)
        #elif action == readlengthsAction:
        #    self.app.read_distributon(row)
        #elif action == contamaction:
        #    self.app.check_contamination()
        elif action == samplereadsaction:
            self.app.get_fasta_reads()
        elif action == normalisefastqaction:
            self.app.normalise_fastq(row)
        elif action == mappingstatsAction:
            self.app.mapping_stats(row)
        #elif action == plotbamAction:
        #    self.app.show_bam_viewer(row)
        elif action == removeAction:
            self.deleteRows(rows)
        elif action == removebamAction:
            self.deleteBamFiles(rows)
        elif action == movefastqAction:
            self.app.movefastq()
        elif action == copyAction:
            self.copy()
        elif action == exportAction:
            self.exportTable()
        #elif action == colorbyAction:
        #    self.colorByColumn()
        return

    def edit(self, index, trigger, event):
        """Override edit to disable editing of columns"""

        if index.column() < 20:
            return False
        else:
            QTableView.edit(self, index, trigger, event)
        return True

    def resizeColumns(self):

        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        #self.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        return

    def deleteRows(self, rows):
        """Remove the sample and it's bam file"""

        answer = QMessageBox.question(self, 'Delete Entry?',
                             'Are you sure? This will remove the alignments and calls but not the fastq file.',
                             QMessageBox.Yes, QMessageBox.No)
        if answer == QMessageBox.No:
            return
        df = self.model.df
        idx = df.index[rows]
        #remove any bam files first
        self.deleteBamFiles(rows, ask=False)
        #self.deleteVariantCalls(rows)
        self.model.df = df.drop(idx)
        self.refresh()
        return

    def deleteVariantCalls(self, rows):
        """Delete folder with var calls for samples"""

        df = self.model.df
        samples = list(df.index[rows])
        #app.delete_variant_calls(samples)
        return

    def deleteBamFiles(self, rows, ask=True):
        """Remove bam files for selected rows if present"""

        if ask == True:
            answer = QMessageBox.question(self, 'Delete Entry?',
                                 'Are you sure? This removes bam files only.',
                                 QMessageBox.Yes, QMessageBox.No)
            if answer == QMessageBox.No:
                return
        df = self.model.df
        if not 'bam_file' in df.columns:
            return
        idx = df.index[rows]
        for file in df.loc[idx].bam_file:
            if pd.isnull(file):
                continue
            if os.path.exists(file):
                os.remove(file)
            index=file+'.bai'
            if os.path.exists(index):
                os.remove(index)
        df.loc[idx,['bam_file','coverage','meandepth']] = np.nan
        self.refresh()
        return

class ResultsTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        self.app = app
        self.setWordWrap(False)

    def addActions(self, event, row, column):

        menu = self.menu
        showSequencesAction = menu.addAction("Show sequences")
        showAlignmentAction = menu.addAction("Alignment for all samples")

        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == showSequencesAction:
            self.app.show_fasta_sequences(row)
        elif action == showAlignmentAction:
            self.app.show_gene_alignment(row)
        return

class Verticalheaderview(QHeaderView):
    """Vertical header view"""
    def __init__(self, parent=None, height=100):
        super().__init__(Qt.Horizontal, parent)
        self._font = QFont(core.FONT)
        self._metrics = QFontMetrics(self._font)
        self._margin = 10
        self.height = height
        return

    def paintSection(self, painter, rect, index):
        # Get the header text
        data = self._get_data(index)
        if not data:  # Skip empty sections
            return

        painter.save()
        # Translate painter to the center of the rect, then rotate
        painter.translate(rect.x() + rect.width() / 2, rect.y() + rect.height())
        painter.rotate(-90)
        # Set font and draw text
        painter.setFont(self._font)
        # Draw text using QRect for proper positioning
        #text_rect = QtCore.QRect(0, 0, rect.height(), rect.width())
        text_rect = QtCore.QRect(-int(rect.height() / 2), -int(rect.width() / 2),
                                 rect.height()*2, rect.width())

        painter.drawText(
            text_rect,
            Qt.AlignCenter,
            data
        )
        painter.restore()  # Restore the painter's state

    def sizeHint(self):
        # Calculate header height based on the tallest text
        #max_height = self._get_text_width() + 2 * 1000
        #return QtCore.QSize(max_height, super().sizeHint().height())
        return QtCore.QSize(20, self.height+10)

    def _get_text_width(self):
        # Calculate the maximum width of all header texts
        return max(self._metrics.width(str(self._get_data(i)))
                   for i in range(self.model().columnCount()))

    def _get_data(self, index):
        # Fetch header data from the model
        return self.model().headerData(index, Qt.Horizontal, Qt.DisplayRole)

class SNPTableModel(DataFrameModel):

    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.colors = {'A':'lightblue', 'T':'#ff704d', 'C':'#5cd65c',
                        'G':'yellow','N':'gray'}
        self.df = dataframe

    def data(self, index, role):
        """see https://www.pythonguis.com/tutorials/qtableview-modelviews-numpy-pandas/"""

        colors = self.colors
        i = index.row()
        j = index.column()
        rowname = self.df.index[i]
        value = self.df.iloc[i, j]
        if role == QtCore.Qt.DisplayRole:
            return value
        elif role == QtCore.Qt.BackgroundRole:
            if (isinstance(value, str)):
                clr = QColor(colors[value])
                if rowname == 'ref':
                    return clr.darker(140)
                else:
                    return clr
        elif role == QtCore.Qt.TextAlignmentRole:
            return QtCore.Qt.AlignCenter

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class SNPTable(DataFrameTable):
    """
    Table for SNP alignment view.
    See docs for format of file.
    """
    def __init__(self, parent=None, app=None, dataframe=None, **kwargs):
        DataFrameTable.__init__(self, parent, dataframe, **kwargs)
        tm = SNPTableModel(dataframe)
        self.setModel(tm)
        self.model = tm
        self.app = app
        h = self.getMaxHeight()
        headerview = Verticalheaderview(parent=self, height=h)
        self.setHorizontalHeader(headerview)
        hh = self.horizontalHeader()
        hh.setDefaultSectionSize(20)
        hh.setVisible(True)
        #hh.customContextMenuRequested.connect(self.columnHeaderMenu)
        self.transposed = False
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        tm = SNPTableModel(df)
        self.setModel(tm)
        self.model = tm
        self.refresh()
        return

    def addActions(self, event, row, column):
        """Right click action menu"""

        df = self.model.df
        menu = self.menu
        #gotoPositionAction = menu.addAction("Find Position")
        showAlignmentAction = menu.addAction("View Alignments")
        #transposeAction = menu.addAction("Transpose")
        orderbyPhyloAction = menu.addAction("Order By Phylogeny")
        loadAction = menu.addAction("Load SNP table")
        zoominAction = menu.addAction("Zoom in")
        zoomoutAction = menu.addAction("Zoom out")
        action = menu.exec_(event.globalPos())

        '''if action == transposeAction:
            df = df.T
            self.setDataFrame(df)
            #self.refresh()
            self.transposed = not self.transposed
        elif action == gotoPositionAction:
            items = [str(i) for i in df.columns]
            pos,ok = QInputDialog.getItem(self, "Go to Position","Pos:", items, 0, False)
            if not ok:
                return
            if self.transposed == True:
                idx = df.index.get_loc(pos)
                self.setScrollPosition(idx,0)
            else:
                col = df.columns.get_loc(int(pos))
                self.setScrollPosition(col,0)
        '''
        if action == orderbyPhyloAction:
            self.orderByPhylogeny()

        elif action == showAlignmentAction:
            #print (row, column)
            pos = df.columns[column]
            sample = df.index[row]
            #print (sample,pos)
            sample_df = self.app.fastq_table.model.df
            bam_file = sample_df.loc[sample].bam_file
            label = f'{sample}:{pos}'
            bm = self.app.launch_bam_viewer(bam_file, label, pos-200)
        elif action == loadAction:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="Text Files(*.csv *.txt);;All Files(*.*)" )
            if not filename:
                return
            df = pd.read_csv(filename, sep=' ', index_col=0).T
            self.setDataFrame(df)
        elif action == zoomoutAction:
            self.zoomOut()
        elif action == zoominAction:
            self.zoomIn()
        #elif action == uniquepositionsAction:
        #    self.app.show_unique_positions(row)
        return

    def columnHeaderMenu(self, pos):
        """Custom column header menu"""

        hheader = self.horizontalHeader()
        idx = hheader.logicalIndexAt(pos)
        column = self.model.df.columns[idx]
        menu = QMenu(self)
        sortAction = menu.addAction("Sort \u2193")
        sortDescAction = menu.addAction("Sort \u2191")
        action = menu.exec_(self.mapToGlobal(pos))
        if action == sortAction:
            self.sort(idx)
        elif action == sortDescAction:
            self.sort(idx, ascending=False)
        return

    def orderByPhylogeny(self):
        """Cluster labels by phylogeny tip order"""

        treefile = self.app.treefile
        if not os.path.exists(treefile):
            print ('no tree file found')
            return
        from Bio import Phylo
        tree = Phylo.read(treefile, "newick")
        #get order of tips and use to order dataframe
        df = self.model.df
        tips = [leaf.name for leaf in tree.get_terminals()]
        #tips.remove('ref')
        if self.transposed == True:
            df = df.reindex(tips, axis=1)
        else:
            df = df.reindex(tips)
        self.setDataFrame(df)
        return

class DistMatrixTableModel(DataFrameModel):

    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.df = dataframe
        self.updateColors()

    def data(self, index, role):
        """see https://www.pythonguis.com/tutorials/qtableview-modelviews-numpy-pandas/"""

        if self.df is None:
            return
        i = index.row()
        j = index.column()
        rowname = self.df.index[i]
        value = self.df.iloc[i, j]
        if role == QtCore.Qt.DisplayRole:
            return str(value)
        elif role == QtCore.Qt.BackgroundRole:
            clr = self.colors.iloc[i, j]
            qc = QColor()
            qc.setRgbF(*clr)
            return qc
        elif role == QtCore.Qt.TextAlignmentRole:
            return QtCore.Qt.AlignCenter

    def updateColors(self):
        """Update color mapping for cells"""

        import matplotlib.colors as mcolors
        df = self.df
        if df is None:
            return

        cmap = plt.cm.Blues
        lut = cmap(df.to_numpy().flatten())
        lut = {i:cmap(i) for i in df.to_numpy().flatten()}
        self.colors = df.applymap(lambda x: lut[x])
        #print (self.colors)
        return

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class DistMatrixTable(DataFrameTable):
    """
    Table for a distance matrix.
    """
    def __init__(self, parent=None, app=None, dataframe=None, **kwargs):
        DataFrameTable.__init__(self, parent, dataframe, **kwargs)
        tm = DistMatrixTableModel(dataframe)
        self.setModel(tm)
        self.model = tm
        self.app = app
        h = self.getMaxHeight()
        headerview = Verticalheaderview(height=h)
        self.setHorizontalHeader(headerview)
        hh = self.horizontalHeader()
        hh.setDefaultSectionSize(30)
        self.transposed = False
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        tm = DistMatrixTableModel(df)
        self.setModel(tm)
        self.model = tm
        tm.updateColors()
        return

    def addActions(self, event, row, column):
        """Right click action menu"""

        df = self.model.df
        menu = self.menu
        orderbyPhyloAction = menu.addAction("Order By Phylogeny")
        loadAction = menu.addAction("Load SNP Matrix")
        zoominAction = menu.addAction("Zoom in")
        zoomoutAction = menu.addAction("Zoom out")
        action = menu.exec_(event.globalPos())

        if action == orderbyPhyloAction:
            self.orderByPhylogeny()
        elif action == loadAction:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="Text Files(*.csv *.txt);;All Files(*.*)" )
            if not filename:
                return
            df = pd.read_csv(filename, sep=' ', index_col=0)
            self.setDataFrame(df)
        elif action == zoomoutAction:
            self.zoomOut()
        elif action == zoominAction:
            self.zoomIn()
        return

    def orderByPhylogeny(self):
        """Cluster labels by phylogeny tip order"""

        #load parent phylo tree
        treefile = self.app.treefile
        from Bio import Phylo
        tree = Phylo.read(treefile, "newick")
        #get order of tips and use to order dataframe
        df = self.model.df
        tips = [leaf.name for leaf in tree.get_terminals()]
        tips.remove('ref')
        #print (tips)
        df = df.reindex(tips)
        df = df.reindex(tips,axis=1)
        #print (df)
        self.setDataFrame(df)
        return

    def refresh(self):
        DataFrameTable.refresh(self)
        #self.model.updateColors()

class CSQTableModel(DataFrameModel):

    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.df = dataframe

    def data(self, index, role):

        i = index.row()
        j = index.column()
        rowname = self.df.index[i]
        value = self.df.iloc[i, j]
        if role == QtCore.Qt.DisplayRole:
            return str(value)
        elif role == QtCore.Qt.BackgroundRole:
            if value == 0 or j <= 3:
                return QColor('white')
            else:
                return QColor('blue')
        elif role == QtCore.Qt.TextAlignmentRole:
            if j>3:
                return QtCore.Qt.AlignCenter

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class CSQTable(DataFrameTable):
    """
    Table for a csq matrix.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        tm = CSQTableModel(dataframe)
        self.setModel(tm)
        self.app = app
        headerview = Verticalheaderview()
        self.setHorizontalHeader(headerview)
        hh = self.horizontalHeader()
        hh.setDefaultSectionSize(18)
        #hh.setSectionResizeMode(QHeaderView.Stretch)
        self.setColumnWidth(0, 80)
        self.setColumnWidth(1, 100)
        self.setColumnWidth(2, 100)
        self.setColumnWidth(3, 80)
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        tm = CSQTableModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def addActions(self, event, row, column):
        """Right click action menu"""

        df = self.model.df
        menu = self.menu
        loadAction = menu.addAction("Load Table")
        zoominAction = menu.addAction("Zoom in")
        zoomoutAction = menu.addAction("Zoom out")
        action = menu.exec_(event.globalPos())

        if action == loadAction:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="Text Files(*.csv *.txt);;All Files(*.*)" )
            if not filename:
                return
            df = pd.read_csv(filename, sep=',')
            self.setDataFrame(df)
        elif action == zoomoutAction:
            self.zoomOut()
        elif action == zoominAction:
            self.zoomIn()
        return

class VCFTableModel(DataFrameModel):

    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.colors = {'A':'lightblue', 'T':'#ff704d', 'C':'#5cd65c',
                        'G':'yellow','N':'gray'}
        self.df = dataframe

    def data(self, index, role):
        """see https://www.pythonguis.com/tutorials/qtableview-modelviews-numpy-pandas/"""

        colors = self.colors
        i = index.row()
        j = index.column()
        rowname = self.df.index[i]
        value = self.df.iloc[i, j]
        if role == QtCore.Qt.DisplayRole:
            return str(value)
        elif role == QtCore.Qt.BackgroundRole:
            if (isinstance(value, str)):
                if value in self.colors:
                    return QColor(self.colors[value])

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class VCFTable(DataFrameTable):
    """
    Table for vcf view.
    See docs for format of file.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        tm = VCFTableModel(dataframe)
        self.setModel(tm)
        hh = self.horizontalHeader()
        hh.setDefaultSectionSize(60)
        self.app = app
        self.transposed = False
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        tm = VCFTableModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def addActions(self, event, row, column):
        """Right click action menu"""

        df = self.model.df
        menu = self.menu
        loadAction = menu.addAction("Load VCF")
        zoominAction = menu.addAction("Zoom in")
        zoomoutAction = menu.addAction("Zoom out")
        exportAction = menu.addAction("Export")
        action = menu.exec_(event.globalPos())

        if action == loadAction:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="Text Files(*.vcf *.vcf.gz);;All Files(*.*)" )
            if not filename:
                return
            df = tools.vcf_to_dataframe(filename).set_index('sample')
            self.setDataFrame(df)
        elif action == zoomoutAction:
            self.zoomOut()
        elif action == zoominAction:
            self.zoomIn()
        elif action == exportAction:
            self.exportTable()