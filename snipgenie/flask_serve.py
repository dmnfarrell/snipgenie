#!/usr/bin/env python

"""
    BTBgenie prototype web app.
    Created Mar 2021
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
#from . import tools, widgets, tables
from flask import Flask, render_template, request
from wtforms import Form, TextField, validators, StringField, SelectField, FloatField
import sqlite3
import folium

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))

DATABASE = '../notebooks/'

webapp = Flask(__name__)

class ControlsForm(Form):
    name = SelectField('name', choices=[])
    path = TextField('path', default='results')
    n = TextField('n', default='2')
    #kinds = [(i,i) for i in plotkinds]
    #kind = SelectField('plot kind', choices=kinds)
    #submit = SubmitField()

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    return db

def help_msg():
    msg = '<a>BTBgenie sample web interface</a><br>'
    #msg += '<a href="%s"> see help page</a>' %wikipage
    return msg

def get_tree():
    """show a tree"""
    
    return

def show_dataframe(df, map):
    """Show points from a dataframe"""

    m = map
    for i,r in df.iterrows():
        popup = get_popup(r)
        cm = folium.CircleMarker(location = [r.LONG, r.LAT],
                               radius = 5,
                               weight = 1,
                               popup = popup,
                               color = 'gray',
                               fill_color = 'blue',
                               fill_opacity = 0.5,
                               fill = True)
        cm.add_to(m)
    return

def get_popup(r):
    """Popup html"""

    html = '<div class="container-fluid"> <div class="row">\
            <h4>%s</h4> <p>seq type: %s</p>'\
            '<p>SB: %s </p> <p>closest: %s</p>'\
            '<p>species: %s </p></div></div>'\
              %(r['name'],r['clade'],r['SB'],r['nearest'],r.species)
    return html

def base_map():

    from folium.plugins import MeasureControl
    coordinate = (53.5, -6.5)
    map = folium.Map(location=coordinate,
                    tiles='cartodbpositron',
                    zoom_start=8,
                     width=700, height=600)
    map.add_child(MeasureControl())
    return map

@webapp.route('/')
def index():
    """main index page"""

    path = request.args.get("path")
    msg = help_msg()

    form = ControlsForm()

    map = base_map()
    df = pd.read_csv('../notebooks/wicklow_test.csv')
    show_dataframe(df, map)
    #map = folium_map._repr_html_()
    map.save('templates/map.html')
    return render_template("index.html",form=form,div='',msg=msg)


def main():
    webapp.run(port=5000, debug=True)

if __name__ == '__main__':
	main()
