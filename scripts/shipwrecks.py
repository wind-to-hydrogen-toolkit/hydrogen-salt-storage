#!/usr/bin/env python
# coding: utf-8

# # Shipwrecks in Irish Waters
#
# <https://isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.5131>

import os
from zipfile import ZipFile

import contextily as cx
import matplotlib.pyplot as plt

from src import read_data as rd

# base data download directory
DATA_DIR = os.path.join("data", "heritage")

URL = (
    "https://gsi.geodata.gov.ie/downloads/Marine/Data/Downloads/Shapefiles/"
    "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"
)

FILE_NAME = "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

data3 = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

data3.head()

data3.shape

data3.crs

ax = data3.to_crs(3857).plot(
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    markersize=15,
    color="chartreuse",
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.title("Shipwrecks in Irish Waters")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()
