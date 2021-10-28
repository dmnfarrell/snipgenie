#!/usr/bin/env python

"""
    Phylogeny viewer.
    Created Feb 2021
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
from . import tools, widgets

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module

class PhyloApp(QMainWindow):
    """Tree viewer with Toytree, wraps TreeViewer widget"""
    def __init__(self, filename=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Tree-viewer")
        self.setGeometry(QtCore.QRect(200, 200, 1000, 600))
        self.setMinimumHeight(150)
        self.main = TreeViewer(self)
        self.main.setFocus()
        self.setCentralWidget(self.main)
        if filename != None:
            self.main.load_tree(filename)
        return

class TreeViewer(QWidget):
    """Phylogeny viewer widget with toytree"""
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
        self.update()
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

    def clear(self):

        self.tree = None
        self.update()
        return

    def update(self):
        """Update the plot"""

        import toytree
        import toyplot

        if self.tree == None:
            self.browser.setHtml('')
            self.tipitems.clear()
            return
        if type(self.tree) is toytree.Multitree.MultiTree:
            self.update_multitree()
            return
        #set colors
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

        dlg = widgets.MultipleInputDialog(self, opts, title='Tree Style', width=300)
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


def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snipgenie gui tool')
    parser.add_argument("-f", "--file", dest="filename",default=[],
                        help="input tree file", metavar="FILE")

    args = vars(parser.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = PhyloApp(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
