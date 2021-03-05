#!/usr/bin/env python

"""
    SQL db viewer.
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
import string
from .qt import *
from . import tools, widgets, tables
from PySide2.QtSql import QSqlDatabase, QSqlQuery, QSqlTableModel

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module

class DBViewer(QMainWindow):
    """DB viewer for MYSQL database"""
    def __init__(self):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("DB-viewer")
        self.setGeometry(QtCore.QRect(200, 200, 1000, 600))
        self.setMinimumHeight(150)
        self.main = QWidget(self)
        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.add_widgets()
        self.create_menu(self)

        return

    def save_data(self):
        """Save layers"""

        data = tools.get_attributes(self)
        data['tree'] = self.tree
        return data

    def load_data(self, data):
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

        self.menubar = self.menuBar()
        self.file_menu = QMenu('File', parent)
        self.file_menu.addAction('Connect to DB', self.load_db)
        #self.file_menu.addAction('Load Test Tree', self.test_tree)

        self.menubar.addMenu(self.file_menu)
        self.tools_menu = QMenu('Tools', parent)
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
        self.splitter.setSizes([300,100])
        self.splitter.setStretchFactor(1,0)
        #layout.addWidget(self.main)

        # Set up the view
        self.view = QTableView()
        self.splitter.addWidget(self.view)

        from PySide2.QtWebEngineWidgets import QWebEngineView
        self.browser = QWebEngineView()
        self.browser.setMinimumSize(200,200)
        self.splitter.addWidget(self.browser)

        #toolswidget = QWidget()
        #self.splitter.addWidget(toolswidget)
        #l = QVBoxLayout(toolswidget)

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
        #self.model.setHeaderData(0, Qt.Horizontal, "animal_id")

        self.model.select()
        self.view.setModel(self.model)
        self.view.resizeColumnsToContents()
        return

    def create_report(self):
        """sample sequencing report"""

        
        return

    def help(self,event=None):
        """Open the online documentation"""

        #import webbrowser
        link='https://github.com/dmnfarrell/btbgenie'
        self.browser.setUrl(link)

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
