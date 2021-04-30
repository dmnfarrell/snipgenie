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
import sqlite3
import folium

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))

DATABASE = '../notebooks/'

webapp = Flask(__name__)

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    return db

def help_msg():
    msg = '<a>BTBgenie sample web interface</a><br>'
    #msg += '<a href="%s"> see help page</a>' %wikipage
    return msg

@webapp.route('/')
def index():
    """main index page"""

    path = request.args.get("path")
    msg = help_msg()

    start_coords = (46.9540700, 142.7360300)
    folium_map = folium.Map(location=start_coords, zoom_start=14)
    #map = folium_map._repr_html_()
    folium_map.save('templates/map.html')
    return render_template("index.html",div='',msg=msg)


def main():
    webapp.run(port=5000, debug=True)

if __name__ == '__main__':
	main()
