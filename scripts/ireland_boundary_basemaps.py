#!/usr/bin/env python
# coding: utf-8

# # Ireland basemaps and boundaries

import os
from zipfile import BadZipFile, ZipFile

import contextily as cx
import matplotlib.pyplot as plt

from src import data as rd

# ## NUTS Island of Ireland boundary
#
# <https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts>

# base data download directory
DATA_DIR = os.path.join("data", "boundaries", "NUTS2021")

URL = (
    "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/"
    "ref-nuts-2021-01m.shp.zip"
)

FILE_NAME = "ref-nuts-2021-01m.shp.zip"

OUT_DIR = os.path.join("data", "basemaps")
os.makedirs(OUT_DIR, exist_ok=True)
# basemap cache directory
cx.set_cache_dir(OUT_DIR)

DATA_DIR_TEMP = os.path.join(DATA_DIR, "temp")

os.makedirs(DATA_DIR_TEMP, exist_ok=True)

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR_TEMP)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

DATA_FILE = os.path.join(DATA_DIR_TEMP, "NUTS_RG_01M_2021_4326_LEVL_1.shp.zip")

ZipFile(DATA_FILE).namelist()

nuts = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

nuts.shape

nuts.head()

nuts.crs

nuts = nuts[nuts["NUTS_ID"].str.contains("IE0|UKN")]

nuts

ie = nuts.dissolve()

ie.bounds

ie.plot(
    color="navajowhite",
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    linewidth=0.4,
)

plt.title("Boundary of the Island of Ireland")
plt.text(-8.75, 51.275, "© EuroGeographics for the administrative boundaries")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## Provinces - OSi National Statutory Boundaries - 2019 - Ungeneralised
#
# <https://data.gov.ie/dataset/provinces-osi-national-statutory-boundaries-2019>

URL = (
    "https://data-osi.opendata.arcgis.com/datasets/"
    "559bc3300384413aa0fe93f0772cb7f1_0.zip?"
    "outSR=%7B%22latestWkid%22%3A2157%2C%22wkid%22%3A2157%7D"
)
DATA_DIR = os.path.join("data", "boundaries")
FILE_NAME = "osi-provinces-ungeneralised-2019.zip"
DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

data = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

data

data.crs

data.shape

data = data.dissolve()

data.bounds

data.plot(
    color="navajowhite",
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    linewidth=0.4,
)

plt.title("Ireland boundary from OSi ungeneralised data")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## Basemaps from xyzservices
#
# <https://xyzservices.readthedocs.io/en/stable/gallery.html>

# bounding box limits with a ~50,000 m buffer
bbox = ie.to_crs(rd.CRS).buffer(5e4).envelope
xmin, ymin, xmax, ymax = bbox.to_crs(3857).total_bounds

ax = ie.plot(
    color="navajowhite",
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    linewidth=0.4,
)
bbox.to_crs(ie.crs).envelope.boundary.plot(ax=ax)
plt.title("Boundary of the Island of Ireland")
plt.text(-8.75, 51.275, "© EuroGeographics for the administrative boundaries")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


def download_basemap(source, zoom):
    """
    Download Contextily basemaps for local use
    """

    out_file = os.path.join(OUT_DIR, f"Ireland.{source['name']}.{zoom}.tif")
    if not os.path.isfile(out_file):
        irl = cx.bounds2raster(
            xmin, ymin, xmax, ymax, path=out_file, zoom=zoom, source=source
        )

    axis = bbox.to_crs(3857).boundary.plot(linewidth=0)
    cx.add_basemap(axis, source=out_file)
    plt.tick_params(labelbottom=False, labelleft=False)
    plt.tight_layout()
    plt.show()


download_basemap(cx.providers.CartoDB.Voyager, 6)

download_basemap(cx.providers.Stamen.Terrain, 6)

download_basemap(cx.providers.USGS.USImagery, 6)
