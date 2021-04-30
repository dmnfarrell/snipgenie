#!/usr/bin/env python

"""
    BTBgenie prototype DB SNP typing viewer.
    Created Jan 2021
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
from . import tools, widgets, tables
from PySide2.QtSql import QSqlDatabase, QSqlQuery, QSqlTableModel

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module

class FoliumWidget(QWidget):
    def __init__(self, parent=None, table=None):
        """Customise this and/or doFrame for your widgets"""

        super(FoliumWidget, self).__init__(parent)

        self.setGeometry(QtCore.QRect(200, 200, 500, 300))
        self.setMinimumHeight(150)
        self.add_widgets()
        #self.create_menu(self)
        self.layers = {}
        self.base_map()
        return

    def create_menu(self, parent):
        """Menu bar"""

        self.menubar = QMenuBar(parent)
        self.file_menu = QMenu('File', parent)
        #self.file_menu.addAction('Import', self.load_map)
        self.menubar.addMenu(self.file_menu)
        return

    def add_widgets(self):
        """Add widgets"""

        layout = self.layout = QVBoxLayout()
        self.main = QWidget()
        vbox = QVBoxLayout(self.main)
        layout.addWidget(self.main)
        from PySide2.QtWebEngineWidgets import QWebEngineView
        self.view = QWebEngineView()
        vbox = QVBoxLayout()
        self.setLayout(vbox)
        vbox.addWidget(self.view)
        self.view.setMinimumHeight(500)
        return

    def import_shapefile(self, filename, name=None, color=None):

        import geopandas as gpd
        ext = os.path.splitext(filename)[1]
        if ext == '.zip':
            filename = 'zip://'+filename
        gdf = gpd.read_file(filename)
        if name == None:
            name = os.path.basename(filename)
        #self.addEntry(name, gdf, color, filename)
        self.show_dataframe(gdf)
        return

    def importFile(self):
        """Import shapefile"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(self,"Import File",
                                                  '.',"shapefile (*.shp);;zip (*.zip);;All files (*.*)",
                                                  options=options)
        if not filename:
            return
        self.import_shapefile(filename)
        return

    def show_dataframe(self, df):
        """Show points from a dataframe"""

        import folium
        m = self.map
        for i,r in df.iterrows():
            popup = self.get_popup(r)
            cm = folium.CircleMarker(location = [r.LONG, r.LAT],
                                   radius = 5,
                                   weight = 1,
                                   popup = popup,
                                   color = 'gray',
                                   fill_color = 'blue',
                                   fill_opacity = 0.5,
                                   fill = True)
            cm.add_to(m)

        data = io.BytesIO()
        m.save(data, close_file=False)
        self.view.setHtml(data.getvalue().decode())
        data.close()
        #self.view.reload()
        return

    def get_popup(self, r):
        """Popup html"""

        html = '<h4>%s</h4> <p>seq type: %s</p>'\
                '<p>SB: %s </p> <p>closest: %s</p>'\
                '<p>species: %s </p>'\
                   %(r['name'],r['clade'],r['SB'],r['nearest'],r.species)
        return html

    def base_map(self):
        """show the base map"""

        import folium
        from folium.plugins import MeasureControl
        coordinate = (53.5, -6.5)
        m = folium.Map(
            tiles='cartodbpositron',
            zoom_start=8,
            location=coordinate
        )
        #folium.TileLayer(tiles='Stamen Terrain',name="Stamen Terrain").add_to(m)
        #folium.LayerControl().add_to(m)
        m.add_child(MeasureControl())

        data = io.BytesIO()
        m.save(data, close_file=False)
        self.view.setHtml(data.getvalue().decode())
        data.close()
        self.map = m
        return

class DBViewer(QMainWindow):
    """Sample app for BTBgenie database/mapping interaction"""
    def __init__(self):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("BTBGenie sample viewer")
        self.setGeometry(QtCore.QRect(200, 200, 1200, 800))
        self.setMinimumHeight(150)
        self.main = QWidget(self)
        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.add_widgets()
        self.create_menu(self)
        self.show_map()
        return

    def create_menu(self, parent):
        """Menu bar"""

        self.menubar = self.menuBar()
        self.file_menu = QMenu('File', parent)
        self.file_menu.addAction('Connect to DB', self.load_db)
        self.file_menu.addAction('Quit', self.quit,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menubar.addMenu(self.file_menu)
        self.tools_menu = QMenu('Tools', parent)
        self.tools_menu.addAction('Map View', self.show_map)
        self.menubar.addMenu(self.tools_menu)
        self.help_menu = QMenu('Help', parent)
        self.menubar.addMenu(self.help_menu)
        self.help_menu.addAction('About', self.help)
        return

    def add_widgets(self):
        """Add widgets"""

        vbox = QVBoxLayout(self.main)
        self.splitter = QSplitter()
        vbox.addWidget(self.splitter)
        self.splitter.setSizes([100,300])
        self.splitter.setStretchFactor(0,1)
        #layout.addWidget(self.main)

        self.left_tabs = QTabWidget()
        self.splitter.addWidget(self.left_tabs)
        # Set up the view
        self.dbview = QTableView()
        self.left_tabs.addTab(self.dbview,'db')

        self.right_tabs = QTabWidget()
        self.right_tabs.setMinimumSize(400,200)
        self.splitter.addWidget(self.right_tabs)

        from PySide2.QtWebEngineWidgets import QWebEngineView
        self.browser = QWebEngineView()
        self.browser.setMinimumSize(300,200)
        #self.splitter.addWidget(self.browser)

        #toolswidget = QWidget()
        #self.splitter.addWidget(toolswidget)
        #l = QVBoxLayout(toolswidget)
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        return

    def load_db(self):
        """connect"""

        con = QSqlDatabase.addDatabase("QSQLITE")
        con.setDatabaseName("mbovis_ireland.sqlite")
        print (con.databaseName())
        print (con.connectionName())
        self.load_table()
        return

    def load_table(self):

        self.model = QSqlTableModel(self)
        self.model.setTable("isolates")
        self.model.setEditStrategy(QSqlTableModel.OnFieldChange)
        self.model.select()
        self.dbview.setModel(self.model)
        self.dbview.resizeColumnsToContents()
        #self.df = self.get_table_data()
        return

    def get_table_data(self):

        #self.model.data
        rows = self.model.rowCount()
        columns = self.model.columnCount()

        for i in range(rows):
            for j in range(columns):
                df.loc[i, j] = str(self.model.item(i, j).text())
        return df

    def get_tabs(self):

        n=[]
        for i in range(self.right_tabs.count()):
            n.append(self.right_tabs.tabText(i))
        return n

    def show_map(self):

        if not hasattr(self, 'mapviewer'):
            self.mapviewer = FoliumWidget()
            self.df = pd.read_csv('notebooks/wicklow_test.csv')
            self.mapviewer.show_dataframe(self.df)
        if not 'map' in self.get_tabs():
            idx = self.right_tabs.addTab(self.mapviewer, 'map')
            self.right_tabs.setCurrentIndex(idx)
        return

    def show_tree(self):

        self.tree_viewer()
        filename = os.path.join(self.outputdir,'RAxML_bipartitions.variants')
        self.treeviewer.load_tree(filename)
        self.treeviewer.update()
        return

    def tree_viewer(self):

        from . import phylo
        if not hasattr(self, 'treeviewer'):
            self.treeviewer = phylo.TreeViewer(self)
        if not 'phylogeny' in self.get_tabs():
            idx = self.right_tabs.addTab(self.treeviewer, 'phylogeny')
            self.right_tabs.setCurrentIndex(idx)
        return

    def help(self,event=None):
        """Open the online documentation"""

        #import webbrowser
        link='https://github.com/dmnfarrell/btbgenie'
        self.browser.setUrl(link)
        return

    def quit(self):
        self.close()
        return

def main():
    "Run the application"

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = DBViewer()
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
