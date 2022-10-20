#!/usr/bin/env python

import os
import sys
from .qt import *
from . import phylo

class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.main = QWidget()
        self.setCentralWidget(self.main)
        self.setGeometry(QtCore.QRect(200, 200, 800, 400))
        l = QHBoxLayout(self.main)
        lbl = QLabel("Tests")
        l.addWidget(lbl)
        self.browser = QWebEngineView()        
        #self.browser.setMinimumHeight(500)
        self.browser.update()
        l.addWidget(self.browser)

        #self.treeview = phylo.TidyTreeViewer(self.main)
        #self.treeview.update()
        #l.addWidget(self.treeview)
        self.show()
        return

    def update(self):
        html = '''
        <html lang="en-UK">
        <body>
        <h2> d3 test</h2>
    	<script src="https://d3js.org/d3.v3.min.js"></script>
        <svg>
          <circle style="fill: #69b3a2" stroke="black" cx=50 cy=50 r=40></circle>
        </svg>

        </body>
        </html>
        '''
        self.browser.setHtml(html)
        return

# creating a pyQt5 application
app = QApplication(sys.argv)
# setting name to the application
app.setApplicationName("javascript Tests")
# creating a main window object
window = MainWindow()

# loop
app.exec_()
