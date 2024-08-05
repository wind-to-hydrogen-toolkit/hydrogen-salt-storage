#!/usr/bin/env python
# coding: utf-8

# # Exploration wells
#
# - <https://data.marine.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2171>
# - <https://data.gov.ie/dataset/exploration-wells-in-the-irish-offshore>

# In[1]:


import os
from zipfile import ZipFile

import contextily as cx
import matplotlib.pyplot as plt

from h2ss import data as rd

# In[10]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0


# In[2]:


# base data download directory
DATA_DIR = os.path.join("data", "exploration-wells")

FILE_NAME = "Exploration_Wells_Irish_Offshore.shapezip.zip"

URL = f"https://atlas.marine.ie/midata/EnergyResourcesExploration/{FILE_NAME}"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[3]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# In[4]:


ZipFile(DATA_FILE).namelist()


# In[5]:


wells = rd.read_shapefile_from_zip(data_path=DATA_FILE)


# In[6]:


wells.shape


# In[7]:


wells.head()


# In[8]:


wells.crs


# In[11]:


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
