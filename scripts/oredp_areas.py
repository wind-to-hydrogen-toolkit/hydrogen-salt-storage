#!/usr/bin/env python
# coding: utf-8

# # OREDP Boundaries

import os
from zipfile import ZipFile

import contextily as cx
import matplotlib.pyplot as plt
import seaborn as sns

from src import data as rd

# ## OREDP Assessment Zone
#
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2212>

# base data download directory
DATA_DIR = os.path.join("data", "oredp-zone")

URL = (
    "https://atlas.marine.ie/midata/EnergyResourcesTidal/"
    "OREDP_Assessment_Zone.shapezip.zip"
)

FILE_NAME = "OREDP_Assessment_Zone.shapezip.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

data = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

data

data.shape

data.columns

data.crs

ax = data.to_crs(3857).plot(
    figsize=(7.5, 7.5),
    column="styleLayer",
    alpha=0.85,
    cmap=sns.color_palette("crest", as_cmap=True),
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.title("OREDP Assessment Zone")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## OREDP Study Area
#
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2214>

URL = (
    "https://atlas.marine.ie/midata/EnergyResourcesTidal/"
    "OREDP_Study_Area.shapezip.zip"
)
FILE_NAME = "OREDP_Study_Area.shapezip.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

data2 = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

data2

data2.shape

data2.columns

data2.crs

ax = data2.to_crs(3857).plot(
    figsize=(9, 9),
    column="localId",
    legend=True,
    legend_kwds={"loc": "lower right"},
    cmap="Accent",
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.title("OREDP Study Area")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()
