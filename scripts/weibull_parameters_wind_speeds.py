#!/usr/bin/env python
# coding: utf-8

# # Weibull parameters of wind speeds - 150 m above ground level
#
# - <https://data.gov.ie/dataset/weibull-parameters-wind-speeds-2001-to-2010-150m-above-ground-level>
# - <https://gis.seai.ie/wind/>

# In[1]:


import os
from zipfile import ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import data as rd
from h2ss import functions as fns

# In[2]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0
plt.rcParams["xtick.minor.size"] = 0
plt.rcParams["ytick.minor.size"] = 0


# In[3]:


# base data download directory
DATA_DIR = os.path.join("data", "weibull-parameters-wind-speeds")

FILE_NAME = "Weibull_150m_params_ITM.zip"

URL = f"https://seaiopendata.blob.core.windows.net/wind/{FILE_NAME}"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[4]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# In[5]:


ZipFile(DATA_FILE).namelist()


# In[6]:


weibull_c = rd.read_shapefile_from_zip(
    data_path=os.path.join(DATA_FILE), endswith="c_ITM.shp"
)


# In[7]:


weibull_k = rd.read_shapefile_from_zip(
    data_path=os.path.join(DATA_FILE), endswith="k_ITM.shp"
)


# In[8]:


weibull_c.crs


# In[9]:


weibull_k.crs


# In[10]:


weibull_c.shape


# In[11]:


weibull_k.shape


# In[12]:


weibull_c.columns


# In[13]:


weibull_k.columns


# In[14]:


weibull_c.head()


# In[15]:


weibull_k.head()


# In[16]:


ds, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))


# In[17]:


xmin, ymin, xmax, ymax = extent.total_bounds


# In[18]:


ax = weibull_c.to_crs(3857).plot(
    column="Value",
    cmap="flare",
    figsize=(6, 6),
    legend=True,
    legend_kwds={"label": "c"},
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


# In[19]:


ax = weibull_k.to_crs(3857).plot(
    column="Value",
    cmap="flare",
    figsize=(6, 6),
    legend=True,
    legend_kwds={"label": "k"},
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


# In[20]:


# wind farms in the area
wind_farms = fns.constraint_wind_farm(
    data_path=os.path.join(
        "data", "wind-farms", "marine-area-consent-wind.zip"
    )
)


# In[21]:


# shape of the halite
shape = rd.halite_shape(dat_xr=ds)


# In[22]:


# land boundary
land = rd.read_shapefile_from_zip(
    data_path=os.path.join(
        "data", "boundaries", "osi-provinces-ungeneralised-2019.zip"
    )
)

land = land.dissolve().to_crs(rd.CRS)


# In[23]:


# crop to wind farm and basin extent
extent_wf = gpd.GeoDataFrame(
    geometry=(
        gpd.GeoDataFrame(geometry=wind_farms.dissolve().envelope)
        .overlay(gpd.GeoDataFrame(geometry=extent), how="union")
        .dissolve()
        .envelope
    )
)
weibull_c = weibull_c.to_crs(rd.CRS).overlay(extent_wf, how="intersection")
weibull_k = weibull_k.to_crs(rd.CRS).overlay(extent_wf, how="intersection")


# In[24]:


# crop land boundary from c and k
weibull_c = weibull_c.overlay(land, how="difference")
weibull_k = weibull_k.overlay(land, how="difference")


# In[25]:


def plot_map(df, label):
    """Plotting helper function"""
    plt.figure(figsize=(10, 10))
    ax1 = plt.axes(projection=ccrs.epsg(rd.CRS))

    # add halite boundary - use buffering to smooth the outline
    shape.buffer(1000).buffer(-1000).boundary.plot(
        ax=ax1, color="black", linewidth=2
    )

    # wind farms
    colours = ["lime", "darkslategrey", "deepskyblue"]
    for index, colour in zip(range(len(wind_farms)), colours):
        wind_farms.iloc[[index]].plot(
            ax=ax1,
            hatch="///",
            facecolor="none",
            edgecolor=colour,
            linewidth=2,
            zorder=2,
        )
    legend_handles = [
        mpatches.Patch(
            facecolor="none",
            hatch="////",
            edgecolor=colours[x],
            label=list(wind_farms["name"])[x],
        )
        for x in range(len(wind_farms))
    ]

    legend_handles.append(
        mpatches.Patch(
            facecolor="none",
            edgecolor="black",
            label="Kish Basin halite",
            linewidth=2,
        )
    )

    df.plot(
        column="Value",
        cmap="flare",
        figsize=(6, 6),
        legend=True,
        ax=ax1,
        zorder=1,
        legend_kwds={"label": label},
    )

    cx.add_basemap(
        ax1, crs=rd.CRS, source=cx.providers.CartoDB.Voyager, zoom=10
    )
    ax1.gridlines(
        draw_labels={"bottom": "x", "left": "y"},
        alpha=0.25,
        color="darkslategrey",
    )
    ax1.add_artist(
        ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
    )
    ax1.legend(handles=legend_handles, loc="upper right")

    plt.title(None)
    plt.tight_layout()
    plt.show()


# In[26]:


plot_map(weibull_c, "c")


# In[27]:


plot_map(weibull_k, "k")


# In[28]:


# areas intersecting with wind farms
weibull_c = weibull_c.overlay(wind_farms, how="intersection")
weibull_k = weibull_k.overlay(wind_farms, how="intersection")


# In[29]:


# compute c and k over wind farms
weibull_c = wind_farms.merge(
    weibull_c.dissolve(by="name", aggfunc={"Value": ["min", "max", "mean"]}),
    on="name",
)
weibull_k = wind_farms.merge(
    weibull_k.dissolve(by="name", aggfunc={"Value": ["min", "max", "mean"]}),
    on="name",
)


# In[30]:


weibull_c[["name", ("Value", "min"), ("Value", "max"), ("Value", "mean")]]


# In[31]:


weibull_k[["name", ("Value", "min"), ("Value", "max"), ("Value", "mean")]]
