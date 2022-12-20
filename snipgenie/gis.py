#!/usr/bin/env python

"""
    snipgenie GIS component.
    Created Jan 2021
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
from .qt import *
import pandas as pd
import numpy as np
import pylab as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from . import widgets
import geopandas as gpd

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
qcolors = ['blue','green','crimson','lightblue','gold','burlywood','blueviolet','chartreuse',
            'cadetblue','coral','cornflowerblue','cornsilk','khaki','orange','pink','chocolate',
            'red','lime','mediumvioletred','navy','teal','darkblue','purple','orange',
            'salmon','brown']

class PlotWidget(FigureCanvas):
    def __init__(self, parent=None, figure=None, dpi=100, hold=False):

        if figure == None:
            figure = Figure()
        super(PlotWidget, self).__init__(figure)
        self.setParent(parent)
        self.figure = Figure(dpi=dpi)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)

class Layer(object):
    def __init__(self, gdf, name, crs=None):
        self.gdf = gdf
        self.name = name
        self.filename = None
        self.crs = crs
        self.lw = 1
        self.color = 'white'
        self.ec = 'black'
        self.alpha = .7
        self.column_color = ''
        self.column_label = ''
        self.colormap = ''
        self.pointsize = 50
        self.labelsize = 10
        return

    def save(self):
        """Save if filename present"""

        if self.filename != None:
            self.gdf.to_file(filename)
        return

class GISViewer(QWidget):
    """Geopandas map plotting widget"""

    def __init__(self, parent=None, table=None):
        """Customise this and/or doFrame for your widgets"""

        super(GISViewer, self).__init__(parent)

        self.parent = parent
        self.tablewidget = table
        self.ID = 'Simple GIS'
        self.layers = {}
        self.createWidgets()
        self.createMenu(self)
        return

    def createMenu(self, parent):
        """Main menu"""

        self.menubar = QMenuBar(parent)
        self.file_menu = QMenu('File', parent)
        self.file_menu.addAction('Import Shapefile', self.importFile)
        self.file_menu.addAction('Load Default Map', self.loadDefault)
        self.file_menu.addAction('Load World Map', self.loadWorldMap)
        self.menubar.addMenu(self.file_menu)
        self.layers_menu = QMenu('Layers', parent)
        self.layers_menu.addAction('Clear', self.clear)
        self.menubar.addMenu(self.layers_menu)
        self.tools_menu = QMenu('Tools', parent)
        self.geom_menu = QMenu('Geometry', self.tools_menu)
        self.geom_menu.addAction('Centroid', lambda: self.apply_geometry('centroid'))
        self.geom_menu.addAction('Boundary', lambda: self.apply_geometry('boundary'))
        self.geom_menu.addAction('Envelope', lambda: self.apply_geometry('envelope'))
        self.geom_menu.addAction('Convex hull', lambda: self.apply_geometry('convex_hull'))
        self.geom_menu.addAction('Buffer', lambda: self.apply_geometry('buffer'))
        self.geom_menu.addAction('Simplify', lambda: self.apply_geometry('simplify'))
        self.geom_menu.addAction('Merge Overlapping', self.mergeOverlap)
        self.tools_menu.addAction(self.geom_menu.menuAction())
        self.transform_menu = QMenu('Transform', self.tools_menu)
        self.transform_menu.addAction('Scale', lambda: self.apply_geometry('scale'))
        self.transform_menu.addAction('Rotate', lambda: self.apply_geometry('rotate'))
        self.transform_menu.addAction('Skew', lambda: self.apply_geometry('skew'))
        self.tools_menu.addAction(self.transform_menu.menuAction())
        self.set_menu = QMenu('Set', self.tools_menu)
        self.set_menu.addAction('Union', self.overlay)
        self.set_menu.addAction('Intersection', lambda: self.overlay('intersection'))
        self.set_menu.addAction('Difference', lambda: self.overlay('difference'))
        self.tools_menu.addAction(self.set_menu.menuAction())
        self.analysis_menu = QMenu('Analysis', self.tools_menu)
        self.tools_menu.addAction(self.analysis_menu.menuAction())
        self.analysis_menu.addAction('Distance Matrix', self.distanceMatrix)
        self.menubar.addMenu(self.tools_menu)
        self.options_menu = QMenu('Options', parent)
        self.subplotsaction = QAction('Multiple Subplots', self.options_menu, checkable=True)
        self.options_menu.addAction(self.subplotsaction)
        self.menubar.addMenu(self.options_menu)
        self.menubar.adjustSize()
        return

    def createWidgets(self):
        """Create widgets"""

        layout = self.layout = QVBoxLayout()
        self.main = QWidget()
        self.vbox = vbox = QVBoxLayout(self.main)
        layout.addWidget(self.main)
        #splitter
        self.splitter = QSplitter(self)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setSizes([50,200])
        self.splitter.setStretchFactor(1,0)
        vbox.addWidget(self.splitter)
        self.setLayout(layout)
        #self.plotview = QWebEngineView()
        pframe = QWidget()
        self.splitter.addWidget(pframe)
        self.addPlotWidget(pframe)
        self.top = QWidget()
        self.splitter.addWidget(self.top)
        #layout.addWidget(self.top)
        hbox = QHBoxLayout(self.top)
        self.tree = QTreeWidget()
        self.tree.setHeaderItem(QTreeWidgetItem(["name","file"]))
        self.tree.setColumnWidth(0, 200)
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self.showTreeMenu)
        #self.tree.itemDoubleClicked.connect(handler)
        #self.tree.itemChanged.connect(self.itemClicked)
        hbox.addWidget(self.tree)
        self.toolbar = self.createToolBar(self)
        hbox.addWidget(self.toolbar)
        return

    def addPlotWidget(self, parent):
        """Create mpl plot canvas and toolbar"""

        layout = QVBoxLayout(parent)
        self.canvas = PlotWidget(parent)
        self.fig = self.canvas.figure
        self.ax = self.canvas.ax
        self.toolbar = NavigationToolbar(self.canvas, parent)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        return

    def createToolBar(self, parent):

        items = {'plot': {'action':self.plot,'file':'plot-map'},
                 'moveup': {'action':lambda: self.moveLayer(1),'file':'arrow-up'},
                 'movedown': {'action':lambda: self.moveLayer(-1),'file':'arrow-down'},
                 'delete': {'action':self.delete,'file':'delete'},
                 }
        toolbar = QToolBar("Toolbar")
        toolbar.setOrientation(QtCore.Qt.Vertical)
        widgets.addToolBarItems(toolbar, self, items)
        return toolbar

    def showTreeMenu(self, pos):
        """Show right cick tree menu"""

        item = self.tree.itemAt( pos )
        menu = QMenu(self.tree)
        editAction = menu.addAction("Edit Table")
        propsAction = menu.addAction("Properties")
        colorAction = menu.addAction("Set Color")
        deleteAction = menu.addAction("Delete")
        setfileAction = menu.addAction("Set File")
        action = menu.exec_(self.tree.mapToGlobal(pos))
        if action == editAction:
            self.edit(item)
        elif action == colorAction:
            self.setColor(item)
        elif action == propsAction:
            self.setProperties(item)
        elif action == deleteAction:
            self.delete(item)
        elif action == setfileAction:
            self.setFile(item)
        return

    def loadData(self, data):
        """Load saved layers"""

        self.updateLayers(data['layers'])
        self.setFigure(data['plot'])
        self.update()
        return

    def setFigure(self, figure):
        """Recreate canvas with given figure"""

        self.canvas.figure = figure
        self.fig = self.canvas.figure
        self.ax = self.canvas.ax
        self.canvas.draw()
        return

    def saveData(self):
        """Save layers"""

        data = {}
        data['layers'] = self.layers
        data['plot'] = self.fig
        print (data)
        return data

    def loadWorldMap(self):
        """Load a world map"""

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        world = world[['continent', 'geometry']]
        self.addEntry('world',world)
        self.plot()
        return

    def loadDefault(self):
        """Load a test map"""

        url = 'https://github.com/dmnfarrell/snipgenie/raw/master/maps/ireland_counties.zip?raw=true'
        self.importShapefile(url, 'ireland counties', 'white')
        return

    def importShapefile(self, filename, name=None, color=None):

        ext = os.path.splitext(filename)[1]
        if ext == '.zip':
            filename = 'zip://'+filename
        gdf = gpd.read_file(filename)
        if name == None:
            name = os.path.basename(filename)
        self.addEntry(name, gdf, color, filename)
        self.plot()
        return

    def importFile(self):
        """Import shapefile"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(self,"Import File",
                                                  '.',"shapefile (*.shp);;zip (*.zip);;All files (*.*)",
                                                  options=options)
        if not filename:
            return
        self.importShapefile(filename)
        return

    def setFile(self, item=None):
        """Attach a file to the layer"""

        if item is None:
            item = self.tree.selectedItems()[0]
        name = item.text(0)
        layer = self.layers[name]
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save to File",
                                                  '.',"shapefile (*.shp);;All files (*.*)",
                                                  options=options)
        ext = os.path.splitext(filename)[1]
        if ext != '.shp':
            filename += '.shp'
        layer.gdf.to_file(filename)
        layer.filename = filename
        item.setText(1, filename)
        return

    def addEntry(self, name, gdf, color=None, filename=None):
        """Add geopandas dataframe entry to tree"""

        if name in self.layers:
            name = name+str(len(self.layers))
        item = QTreeWidgetItem(self.tree)
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        item.setText(0, name)
        item.setText(1, filename)
        #add to layers
        new = Layer(gdf, name)
        self.layers[name] = new
        new.filename = filename
        i = len(self.layers)
        if color == None:
            color = qcolors[i]
        new.color = color
        item.setBackground(0 , QBrush(QColor(color)))
        return

    def updateLayers(self, layers):
        """Reload from a dict of layers"""

        for l in layers:
            lyr=layers[l]
            self.addEntry(lyr.name, lyr.gdf, lyr.color)
        return

    def setColor(self, item):

        qcolor = QColorDialog.getColor()
        item.setBackground(0 , qcolor)
        name = item.text(0)
        self.layers[name].color = qcolor.name()
        self.replot()
        return

    def plot(self, evt=None, ax=None, limits={}):
        """Plot maps"""

        subplots = self.subplotsaction.isChecked()
        order = self.getLayerOrder()
        checked = self.getChecked()
        #get the plot frame from parent table widget
        #pf = self.tablewidget.pf
        column = None
        i=1
        self.fig.clear()
        if ax == None:
            if subplots == 1:
                size = len(order)
                nrows = int(round(np.sqrt(size),0))
                ncols = int(np.ceil(size/nrows))
            else:
                ax = self.fig.add_subplot(111)

        for name in order:
            if name not in checked:
                continue
            if subplots == 1:
                ax = self.fig.add_subplot(nrows,ncols,i,label=name)
                ax.set_title(name)
            layer = self.layers[name]
            clr = layer.color
            df = layer.gdf
            column = layer.column_color
            cmap = layer.colormap
            if layer.column_color != '':
                column=layer.column_color
                df.plot(ax=ax,column=column,cmap=cmap,
                        ec=layer.ec,lw=layer.lw,alpha=layer.alpha)
            else:
                df.plot(ax=ax,color=clr,ec=layer.ec,lw=layer.lw,
                    markersize=layer.pointsize,alpha=layer.alpha)
            #labels
            if layer.column_label != '':
                col = layer.column_label
                df.apply(lambda x: ax.annotate(text=x[col],
                    xy=x.geometry.centroid.coords[0], ha='right', fontsize=layer.labelsize),axis=1)
            if name in limits:
                lims = limits[name]
                ax.set_xlim(lims[0])
                ax.set_ylim(lims[1])
            i+=1
            ax.id = name
        plt.tight_layout()
        self.canvas.draw()
        return

    def plot_folium_points(self, gdf):
        """Plot points with folium"""

        m = self.m
        gdf = gdf.set_crs('EPSG:29902').to_crs('EPSG:4632')
        #colors = tools.random_colors(n=len(labels),seed=20)

        print (gdf)
        for i,r in gdf.iterrows():
            x=r.geometry.x
            y=r.geometry.y
            w=0.005
            pts = ((y-w/1.5,x-w),(y+w/1.5,x+w))
            folium.Circle(location=(y,x), radius=400,
                          color=False,fill=True,fill_opacity=0.6,
                          fill_color='blue',tooltip=r.label).add_to(m)


    def plot_folium(self):
        """Plot with folium"""

        order = self.getLayerOrder()
        checked = self.getChecked()
        column = None
        i=1

        fig = Figure(width=900, height=900)
        m = folium.Map(location=[54.1, -7.0], crs='EPSG3857',tiles='Stamen Terrain',
                              width=1250, height=900)#, min_zoom=12)
        style1 = {'fillColor': 'blue', 'color': 'black','weight':2}

        data = io.BytesIO()
        m.save(data, close_file=False)
        self.plotview.setHtml(data.getvalue().decode())
        self.m = m
        layer = self.layers['centroids.shp']
        #self.plot_points(layer.gdf)
        #layer = self.layers['lpis_merged.shp']
        #print (layer.gdf.crs)
        #p = folium.GeoJson(layer.gdf.to_crs('EPSG:3857'),style_function=lambda x:style1)
        #print (p)
        #m.add_child(p)
        return

    def plot_leaflet(self):
        """ipyleaflet plot"""

        from ipyleaflet import Map, basemaps, basemap_to_tiles

        m = Map(basemap=basemaps.CartoDB.Positron, center=(54.1, -7.0), zoom=9, height=800)

        m.save('temp.html')
        with open('temp.html', 'r') as f:
            html = f.read()
            self.plotview.setHtml(html)
        return

    def getPlotLimits(self):

        axes = self.fig.axes
        limits = {}
        for ax in axes:
            limits[ax.id] = (ax.get_xlim(),ax.get_ylim())
        return limits

    def replot(self):
        """Plot after edits to layers"""

        lims = self.getPlotLimits()
        self.plot(limits=lims)
        return

    def moveLayer(self, n=1):
        """Move layer up in tree"""

        l = self.tree.topLevelItemCount()
        item = self.tree.selectedItems()[0]
        row = self.tree.selectedIndexes()[0].row()
        name = item.text(0)
        if row-n >= l or row-n<0:
            return
        self.tree.takeTopLevelItem(row)
        self.tree.insertTopLevelItem(row - n, item)
        self.tree.setCurrentItem(item)
        return

    def getLayerOrder(self):

        order = []
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            order.append(item.text(0))
        return order[::-1]

    def getChecked(self):

        names=[]
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            if item.checkState(0) == QtCore.Qt.CheckState.Checked:
                names.append(item.text(0))
        return names

    def clear(self):
        """Clear all layers"""

        reply = QMessageBox.question(self, 'Clear All',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.No:
            return False
        self.layers = {}
        self.tree.clear()
        return

    def delete(self):
        """Remove layer"""

        item = self.tree.selectedItems()[0]
        row = self.tree.selectedIndexes()[0].row()
        name = item.text(0)
        del self.layers[name]
        self.tree.takeTopLevelItem(row)
        return

    def edit(self, item):
        """edit dataframe in main table"""

        name = item.text(0)
        layer = self.layers[name]

        from . import tables
        table = tables.DataFrameTable(self)
        table.model.df = layer.gdf
        table.refresh()
        self.vbox.addWidget(table)
        return

    def setProperties(self, item):
        """Set selected item properties"""

        colormaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        name = item.text(0)
        layer = self.layers[name]
        crs_vals = ['','EPSG:29990']
        cols = ['']+list(layer.gdf.columns)
        cols.remove('geometry')
        clrs = ['black','white']+qcolors
        opts = {'name':{'type':'entry','default':layer.name},
                'crs': {'type':'combobox','default':layer.crs,'items':crs_vals,'label':'CRS'},
                'line width': {'type':'spinbox','default':layer.lw,'range':(0,15)},
                'edge color': {'type':'combobox','default':layer.ec,'items':clrs},
                'pointsize': {'type':'spinbox','default':layer.pointsize,'range':(50,500)},
                'alpha': {'type':'spinbox','default':layer.alpha,'range':(0.1,1),'interval':0.1},
                'colorby': {'type':'combobox','default':layer.column_color,'items':cols,
                'label':'color by'},
                'labelby': {'type':'combobox','default':layer.column_label,'items':cols,
                'label':'labels'},
                'labelsize':  {'type':'spinbox','default':layer.labelsize,'range':(5,40)},
                'cmap': {'type':'combobox','default':layer.colormap,'items':colormaps},
                }
        dlg = widgets.MultipleInputDialog(self, opts, title='Layer Properties', width=200)
        dlg.exec_()
        if not dlg.accepted:
            return False

        layer.name = dlg.values['name']
        layer.lw = dlg.values['line width']
        layer.ec = dlg.values['edge color']
        layer.pointsize = dlg.values['pointsize']
        layer.alpha = dlg.values['alpha']
        layer.column_color = dlg.values['colorby']
        layer.column_label = dlg.values['labelby']
        layer.colormap = dlg.values['cmap']
        layer.labelsize = dlg.values['labelsize']
        item.setText(0, name)
        self.replot()
        return

    def export(self, item):

        return

    def overlay(self, how='union'):
        """Find difference"""

        items = self.tree.selectedItems()[:2]
        print (items)
        names = [i.text(0) for i in items]
        df1 = self.layers[names[0]].gdf
        df2 = self.layers[names[1]].gdf
        new = gpd.overlay(df1, df2, how=how)
        self.addEntry(names[0]+'_'+names[1]+'_'+how, new)
        self.replot()
        return

    def apply_geometry(self, func):
        """Apply function to geoseries"""

        params = {'buffer':
                  {'distance':{'type':'spinbox','default':0.1,'range':(1,100),'interval':0.1}},
                  'simplify':
                  {'tolerance':{'type':'spinbox','default':0.1,'range':(1,10),'interval':0.1}},
                  'scale':
                  {'xfact':{'type':'spinbox','default':0.1,'range':(1,10),'interval':0.1},
                   'yfact':{'type':'spinbox','default':0.1,'range':(1,10),'interval':0.1}},
                  'rotate':
                  {'angle':{'type':'spinbox','default':0.1,'range':(1,180),'interval':0.2}},
                  'skew':
                  {'xs':{'type':'spinbox','default':0.1,'range':(1,180),'interval':0.2},
                   'ys':{'type':'spinbox','default':0.1,'range':(1,180),'interval':0.2}, },
                 }
        if func in params:
            opts = params[func]
            dlg = widgets.MultipleInputDialog(self, opts, title=func, width=200)
            dlg.exec_()
            if not dlg.accepted:
                return False

        item = self.tree.selectedItems()[0]
        name = item.text(0)
        layer = self.layers[name]
        if func in params:
            new = getattr(layer.gdf.geometry, func)(**dlg.values)
        else:
            new = getattr(layer.gdf.geometry, func)
        new = gpd.GeoDataFrame(geometry=new)
        self.addEntry(name+'_%s' %func, new)
        self.replot()
        return

    def mergeOverlap(self):

        item = self.tree.selectedItems()[0]
        name = item.text(0)
        layer = self.layers[name]
        new = merge_overlap(layer.gdf)
        self.addEntry(name+'_merged', new)
        return

    def distanceMatrix(self):
        """Get dist matrix"""

        item = self.tree.selectedItems()[0]
        name = item.text(0)
        layer = self.layers[name]

        cols = list(layer.gdf.columns)
        opts = {'index':{'type':'combobox','default':cols[0],'items':cols}}
        dlg = widgets.MultipleInputDialog(self, opts, title='distance matrix', width=200)
        dlg.exec_()
        if not dlg.accepted:
            return False

        X = distance_matrix(layer.gdf, index=dlg.values['index'])
        table = self.tablewidget.table
        table.model.df = X
        table.refresh()
        return
