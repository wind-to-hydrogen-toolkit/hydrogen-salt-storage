#!/usr/bin/env python
# coding: utf-8

# # Kish Basin halite data
#
# <https://hyss.ie/>

# In[ ]:


import glob
import os
from textwrap import wrap
from zipfile import BadZipFile, ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from geocube.api.core import make_geocube
from shapely.geometry import Polygon

from h2ss import data as rd

# In[ ]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0


# In[ ]:


# base data download directory
DATA_DIR = os.path.join("data", "kish-basin")

FILE_NAME = "Kish-Basin-dat-files.zip"

URL = f"https://hyss.ie/wp-content/uploads/2023/07/{FILE_NAME}"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[ ]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# In[ ]:


ZipFile(DATA_FILE).namelist()


# In[ ]:


# extract archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)


# ## Map extent

# In[ ]:


with open(
    os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv"),
    encoding="utf-8",
) as f:
    print(f.read())


# In[ ]:


CRS = 23029


# In[ ]:


# create extent polygon
extent = pd.read_csv(
    os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv"), skiprows=2
)
extent = gpd.GeoSeries(
    Polygon(
        [
            (extent[" X"][0], extent[" Y"][0]),
            (extent[" X"][1], extent[" Y"][1]),
            (extent[" X"][2], extent[" Y"][2]),
            (extent[" X"][3], extent[" Y"][3]),
        ]
    ),
    crs=CRS,
)


# In[ ]:


extent.bounds


# In[ ]:


extent.crs


# In[ ]:


plt.figure(figsize=(7, 7))
ax = plt.axes(projection=ccrs.epsg(3857))
plt.xlim(-1.2e6, -0.55e6)
plt.ylim(6.65e6, 7.475e6)
extent.to_crs(3857).boundary.plot(ax=ax, color="darkslategrey")
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)
plt.title("Kish GIS Map Extent")
plt.tight_layout()
plt.show()


# ## XYZ data

# In[ ]:


def read_dat_file(dat_path: str, dat_crs):
    """
    Read XYZ data layers into an Xarray dataset
    """

    gdf = {}
    for dat_file in glob.glob(os.path.join(dat_path, "*.dat")):
        # read each layer as individual dataframes into a dictionary
        gdf[os.path.split(dat_file)[1][:-4]] = pd.read_fwf(
            dat_file, header=None, names=["X", "Y", "Z"]
        )

        # assign layer name to a column
        gdf[os.path.split(dat_file)[1][:-4]]["data"] = os.path.split(dat_file)[
            1
        ][:-4]
    dat_file = glob.glob(os.path.join(dat_path, "*.dat"))[-1]

    # find data resolution
    gdf_xr = (
        gdf[os.path.split(dat_file)[1][:-4]].set_index(["X", "Y"]).to_xarray()
    )
    resx = gdf_xr["X"][1] - gdf_xr["X"][0]
    resy = gdf_xr["Y"][1] - gdf_xr["Y"][0]

    # combine dataframes
    gdf = pd.concat(gdf.values())

    # convert dataframe to geodataframe
    gdf["wkt"] = (
        "POINT (" + gdf["X"].astype(str) + " " + gdf["Y"].astype(str) + ")"
    )
    gdf = gpd.GeoDataFrame(
        gdf, geometry=gpd.GeoSeries.from_wkt(gdf["wkt"]), crs=dat_crs
    )
    gdf.drop(columns=["wkt", "X", "Y"], inplace=True)

    # convert to Xarray dataset
    gdf = make_geocube(
        vector_data=gdf,
        resolution=(-abs(resy), abs(resx)),
        align=(abs(resy / 2), abs(resx / 2)),
        group_by="data",
    )

    return gdf


# In[ ]:


ds = read_dat_file(DATA_DIR, dat_crs=CRS)


# In[ ]:


ds


# In[ ]:


ds.rio.crs


# In[ ]:


ds.rio.resolution()


# In[ ]:


ds.rio.bounds()


# In[ ]:


xmin, ymin, xmax, ymax = extent.total_bounds


# In[ ]:


def plot_maps(plot_data):
    """
    Create facet contour plots of the Xarray dataset
    """

    fig = plot_data["Z"].plot.contourf(
        col="data",
        cmap="jet",
        col_wrap=2,
        robust=True,
        levels=15,
        subplot_kws={"projection": ccrs.epsg(CRS)},
        xlim=(xmin, xmax),
        ylim=(ymin, ymax),
        cbar_kwargs={"aspect": 20, "pad": 0.02},
    )

    # add a basemap
    basemap = cx.providers.CartoDB.PositronNoLabels
    for n, axis in enumerate(fig.axs.flat):
        cx.add_basemap(axis, crs=CRS, source=basemap, attribution=False)
        # add attribution for basemap tiles
        if n == 2:
            axis.text(xmin, ymin - 2500, basemap["attribution"], fontsize=8)

    # assign titles this way to prevent truncation/overflow
    for axis, title in zip(fig.axs.flat, plot_data["data"].values):
        axis.set_title("\n".join(wrap(title, 34)), fontsize=10)
    plt.show()


# ### Halite thickness

# In[ ]:


plot_maps(ds.sel(data=[x for x in ds["data"].values if "Thickness XYZ" in x]))


# ### Halite thickness - zones of interest

# In[ ]:


plot_maps(ds.sel(data=[x for x in ds["data"].values if "Zone" in x]))


# ### Halite base depth

# In[ ]:


plot_maps(ds.sel(data=[x for x in ds["data"].values if "Base" in x]))


# ### Halite top depth

# In[ ]:


plot_maps(ds.sel(data=[x for x in ds["data"].values if "Top Depth" in x]))


# ### Halite top TWT (two-way thickness)

# In[ ]:


plot_maps(ds.sel(data=[x for x in ds["data"].values if "Millisecond" in x]))
