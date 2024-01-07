#!/usr/bin/env python
# coding: utf-8

# # Frequently used shipping routes (300 gross tonnes and above)
#
# <https://data.gov.ie/dataset/frequently-used-routes-300-gross-tonnes-and-above1>

import os
from zipfile import ZipFile

import contextily as cx
import matplotlib.pyplot as plt

from hydrogen_salt_storage import data as rd

# base data download directory
DATA_DIR = os.path.join("data", "shipping")

URL = (
    "https://data-housinggovie.opendata.arcgis.com/datasets/"
    "housinggovie::frequently-used-routes-300-gross-tonnes-and-above.zip?"
    "outSR=%7B%22latestWkid%22%3A3857%2C%22wkid%22%3A102100%7D"
)

FILE_NAME = "shipping_frequently_used_routes.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

data = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

data.shape

data.head()

data.crs

data = data.sort_values("weight", ascending=False)

data.weight.unique()

data.key_.unique()

data["weight"] = data["weight"].astype(str)

data["legend"] = data.weight + ": " + data.key_

ax = data.plot(
    figsize=(7.5, 7.5),
    column="legend",
    legend=True,
    cmap="tab20",
    legend_kwds={"loc": "upper left"},
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)

plt.title("Frequently Used Routes (300 gross tonnes and above)")

plt.tick_params(labelbottom=False, labelleft=False)
# ax.legend(list(data.key_.unique()))
plt.tight_layout()
plt.show()
