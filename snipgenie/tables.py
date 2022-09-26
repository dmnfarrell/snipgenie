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

class ColumnHeader(QHeaderView):
    def __init__(self):
        super(QHeaderView, self).__init__()
        return

class DataFrameWidget(QWidget):
    """Widget containing a tableview and statusbar"""
    def __init__(self, parent=None, table=None, toolbar=False):

        super(DataFrameWidget, self).__init__()
        l = self.layout = QGridLayout()
        l.setSpacing(2)
        self.setLayout(self.layout)
        #self.table = DataFrameTable(self, dataframe)
        self.table = table
        l.addWidget(self.table, 1, 1)
        if toolbar==True:
            self.createToolbar()
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        self.updateStatusBar()
        return

    def createToolbar(self):
        return

    def updateStatusBar(self):
        """Update the table details in the status bar"""

        if not hasattr(self, 'size_label'):
            return
        df = self.table.model.df
        meminfo = self.table.getMemory()
        s = '{r} rows x {c} columns | {m}'.format(r=len(df), c=len(df.columns),m=meminfo)
        self.size_label.setText(s)
        return

class DataFrameTable(QTableView):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, dataframe=None, plotter=None, fontsize=10, *args):

        QTableView.__init__(self)
        self.clicked.connect(self.showSelection)
        #self.doubleClicked.connect(self.handleDoubleClick)
        self.setSelectionBehavior(QTableView.SelectRows)
        #self.setSelectionBehavior(QTableView.SelectColumns)
        #self.horizontalHeader = ColumnHeader()
        header = self.horizontalHeader()
        vh = self.verticalHeader()
        vh.setVisible(True)
        vh.setDefaultSectionSize(28)
        hh = self.horizontalHeader()
        hh.setVisible(True)
        hh.setSectionsMovable(True)
        hh.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        hh.customContextMenuRequested.connect(self.columnHeaderMenu)
        hh.setSelectionBehavior(QTableView.SelectColumns)
        #hh.sectionClicked.connect(self.columnClicked)

        self.setDragEnabled(True)
        self.viewport().setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.resizeColumnsToContents()
        self.setCornerButtonEnabled(True)

        self.font = QFont("Arial", fontsize)
        #print (fontsize)
        self.setFont(self.font)
        tm = DataFrameModel(dataframe)
        self.setModel(tm)
        self.model = tm
        self.setWordWrap(False)
        self.setCornerButtonEnabled(True)
        self.plotview = plotter
        return

    def createToolbar(self):

        self.toolbar = ToolBar(self)
        self.layout.addWidget(self.toolbar, 1, 2)

    def setDataFrame(self, df):

        tm = DataFrameModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def getDataFrame(self):
        return self.model.df

    def load(self, filename=None):
        return

    def save(self):
        return

    def copy(self):
        return

    def paste(self):
        return

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

    def handleDoubleClick(self, item):

        cellContent = item.data()
        if item.column() != 0:
            return
        return

    def columnClicked(self, col):

        hheader = self.horizontalHeader()
        df = self.model.df
        self.model.df = df.sort_values(df.columns[col])
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
        #model = self.model
        menu = QMenu(self)
        sortAction = menu.addAction("Sort \u2193")
        sortDescAction = menu.addAction("Sort \u2191")
        deleteColumnAction = menu.addAction("Delete Column")
        renameColumnAction = menu.addAction("Rename Column")
        addColumnAction = menu.addAction("Add Column")
        plotAction = menu.addAction("Histogram")

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
        self.addActions(event, row)
        return

    def addActions(self, event, row):

        menu = self.menu
        copyAction = menu.addAction("Copy")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == copyAction:
            self.copy()
        elif action == exportAction:
            self.exportTable()
        return

    def setIndex(self):
        return

    def copy(self):

        self.model.df
        return

    def refresh(self):

        self.model.beginResetModel()
        self.model.dataChanged.emit(0,0)
        self.model.endResetModel()

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

    def plotHist(self, column):

        df = self.model.df
        d = df[column]
        #fig,ax = plt.subplots(1,1, figsize=(7,5), dpi=120)
        d.hist(bins=12, ax=self.plotview.ax)
        self.plotview.redraw()
        return

class DataFrameModel(QtCore.QAbstractTableModel):
    def __init__(self, dataframe=None, *args):
        super(DataFrameModel, self).__init__()
        if dataframe is None:
            self.df = pd.DataFrame()
        else:
            self.df = dataframe
        self.bg = '#F4F4F3'
        return

    def update(self, df):
        #print('Updating Model')
        self.df = df

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.columns.values)

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
                elif type(value) == np.float:
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
            return QColor(self.bg)

    def headerData(self, col, orientation, role=QtCore.Qt.DisplayRole):

        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return str(self.df.columns[col])
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            return self.df.index[col]
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
        #if role == QtCore.Qt.DisplayRole:
        #    return value
        color = QColor(self.bg)
        if role == QtCore.Qt.BackgroundRole:
            #color warnings
            if colname == 'meandepth':
                if value < 20:
                    color = QColor('#ff704d')
            elif colname == 'coverage':
                if value < 98:
                    color = QColor('#ff704d')
            else:
                color = QColor(self.bg)
            return color
        else:
            return super(SampleTableModel, self).data(index, role)

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class SampleTable(DataFrameTable):
    """
    QTableView for files/samples view.
    """
    def __init__(self, parent=None, app=None, dataframe=None, plotter=None, *args):
        DataFrameTable.__init__(self)
        self.app = app
        self.setWordWrap(False)
        header = self.horizontalHeader()
        tm = SampleTableModel(dataframe)
        self.setModel(tm)
        self.plotview = plotter
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        tm = SampleTableModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def addActions(self, event, row):

        menu = self.menu
        detailsAction = menu.addAction("Sample Details")
        fastqqualityAction = menu.addAction("Quality Summary")
        readlengthsAction = menu.addAction("Read Length Distribution")
        plotbamAction = menu.addAction("Show Read Alignments")
        mappingstatsAction = menu.addAction("Mapping Statistics")
        contamaction = menu.addAction('Check Contamination')
        fastareadsaction = menu.addAction('Sample Sequences')
        removeAction = menu.addAction("Remove Selected")
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
        elif action == readlengthsAction:
            self.app.read_distributon(row)
        elif action == contamaction:
            self.app.check_contamination()
        elif action == fastareadsaction:
            self.app.get_fasta_reads()
        elif action == mappingstatsAction:
            self.app.mapping_stats(row)
        elif action == plotbamAction:
            #print (row)
            self.app.show_bam_viewer(row)
        elif action == removeAction:
            self.deleteRows(rows)
        elif action == exportAction:
            self.exportTable()
        return

    def edit(self, index, trigger, event):
        """Override edit to disable editing of columns"""

        if index.column() < 20:
            return False
        else:
            QTableView.edit(self, index, trigger, event)
        return True

    def refresh(self):
        DataFrameTable.refresh(self)

    def resizeColumns(self):

        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        #self.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        return

    def deleteRows(self, rows):

        answer = QMessageBox.question(self, 'Delete Entry?',
                             'Are you sure? This will not remove the sample file.',
                             QMessageBox.Yes, QMessageBox.No)
        if answer == QMessageBox.No:
            return
        idx = self.model.df.index[rows]
        self.model.df = self.model.df.drop(idx)
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

    def addActions(self, event, row):

        menu = self.menu
        showSequencesAction = menu.addAction("Show sequences")
        showAlignmentAction = menu.addAction("Alignment for all samples")

        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == showSequencesAction:
            self.app.show_fasta_sequences(row)
        elif action == showAlignmentAction:
            self.app.show_gene_alignment(row)
        return

class MyHeaderView(QHeaderView):

    def __init__(self, parent=None):
        super().__init__(Qt.Horizontal, parent)
        self._font = QFont("Arial", 12)
        self._metrics = QFontMetrics(self._font)
        self._descent = self._metrics.descent()
        self._margin = 10

    def paintSection(self, painter, rect, index):
        data = self._get_data(index)
        painter.rotate(-90)
        painter.setFont(self._font)
        painter.drawText(- rect.height() + self._margin,
                         rect.left() + (rect.width() + self._descent) / 2, data)

    def sizeHint(self):
        return QtCore.QSize(0, self._get_text_width() + 2 * self._margin)

    def _get_text_width(self):
        return max([self._metrics.width(self._get_data(i))
                    for i in range(0, self.model().columnCount())])

    def _get_data(self, index):
        return self.model().headerData(index, self.orientation())

class SNPTableModel(DataFrameModel):

    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.colors = {'A':'lightblue', 'T':'#ff704d', 'C':'#5cd65c', 'G':'yellow'}
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
                if rowname == 'ref':
                    return QColor('gray')
                else:
                    return QColor(colors[value])

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled

class SNPTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, app=None, dataframe=None, *args):
        DataFrameTable.__init__(self, parent, dataframe)
        tm = SNPTableModel(dataframe)
        self.setModel(tm)
        self.app = app
        headerview = MyHeaderView()
        self.setHorizontalHeader(headerview)
        hh = self.horizontalHeader()
        hh.setDefaultSectionSize(18)
        return

    def setDataFrame(self, df):
        """Override to use right model"""

        tm = SNPTableModel(df)
        self.setModel(tm)
        self.model = tm
        return

    def addActions(self, event, row):
        menu = self.menu
        uniquepositionsAction = menu.addAction("View unique positions")

        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == uniquepositionsAction:
            self.app.show_unique_positions(row)

        return
