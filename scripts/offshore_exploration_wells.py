#!/usr/bin/env python
# coding: utf-8

# # Exploration wells
#
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2171>

import os
from zipfile import ZipFile

import contextily as cx
import matplotlib.pyplot as plt

from src import data as rd

# base data download directory
DATA_DIR = os.path.join("data", "exploration-wells")

URL = (
    "https://atlas.marine.ie/midata/EnergyResourcesExploration/"
    "Exploration_Wells_Irish_Offshore.shapezip.zip"
)

FILE_NAME = "Exploration_Wells_Irish_Offshore.shapezip.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

wells = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

wells.shape

wells.head()

wells.crs

ax = wells.to_crs(3857).plot(
    column="AREA",
    legend=True,
    cmap="tab20b",
    figsize=(7.5, 7.5),
    legend_kwds={"loc": "upper right"},
    linewidth=0.5,
    edgecolor="darkslategrey",
)
plt.xlim(-1.6e6, -0.2e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=6)

plt.title("Exploration Wells in the Irish Offshore (1970-2019)")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()
