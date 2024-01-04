#!/usr/bin/env python
# coding: utf-8

# # Bathymetry
#
# <https://isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.858>

import os
from zipfile import BadZipFile, ZipFile

import contextily as cx
import matplotlib.pyplot as plt
import rioxarray as rxr

from src import read_data as rd

# base data download directory
DATA_DIR = os.path.join("data", "bathymetry")

FILE_NAME = "IE_GSI_MI_Bathymetry_25m_IE_Waters_WGS84_LAT_TIFF.zip"

URL = (
    "https://gsi.geodata.gov.ie/downloads/Marine/Data/Downloads/"
    "LatestEntireAreaMerge/" + FILE_NAME
)

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

# # extract the archive
# try:
#     z = ZipFile(DATA_FILE)
#     z.extractall(DATA_DIR)
# except BadZipFile:
#     print("There were issues with the file", DATA_FILE)

data = rxr.open_rasterio(DATA_FILE.split(".")[0] + ".tif", chunks="auto")

data

data.rio.crs

data.rio.bounds()

data.rio.resolution()

# read Kish Basin extent
_, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))

extent.crs

# convert the extent's CRS to match the raster's
extent_ = extent.to_crs(data.rio.crs)

extent_.bounds

# data.rio.clip(extent_)

# clip to extent
data = rxr.open_rasterio(
    DATA_FILE.split(".")[0] + ".tif", chunks="auto", masked=True
).rio.clip()
# data = rxr.open_rasterio(FILE_PATH, cache=False, masked=True).rio.clip(
#     ie["geometry"], from_disk=True
# )
