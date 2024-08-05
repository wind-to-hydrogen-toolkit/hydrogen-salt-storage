#!/usr/bin/env python
# coding: utf-8

# # KIS-ORCA subsea cables
#
# <https://kis-orca.org/>
#
# **IMPORTANT:** There may be some incorrect name assignments as a simple
# backfill and multiline conversion was used and no further checks were done

# In[1]:


import os
import subprocess
from zipfile import BadZipFile, ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

from h2ss import data as rd

# In[2]:


# base data download directory
DATA_DIR = os.path.join("data", "subsea-cables")

FILE_NAME = "Olex_KIS-ORCA-v2023.zip"

URL = f"https://kis-orca.org/wp-content/uploads/2020/12/{FILE_NAME}"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[3]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0


# In[4]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# ## Read data

# In[5]:


# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)


# In[6]:


DATA_FILE = os.path.join(DATA_DIR, ZipFile(DATA_FILE).namelist()[0])


# In[7]:


# extract gz file using gzip
# https://www.gnu.org/software/gzip/
subprocess.run(["gzip", "-d", "<", DATA_FILE, ">", DATA_FILE[:-3]])


# In[8]:


DATA_FILE = DATA_FILE[:-3]


# In[9]:


for n, line in enumerate(open(DATA_FILE, "r", encoding="ISO-8859-15")):
    if n < 11:
        print(line[:-1])


# In[10]:


# coordinates
with open(f"{DATA_FILE}_data.txt", "w", encoding="utf-8") as outfile:
    for n, line in enumerate(open(DATA_FILE, "r", encoding="ISO-8859-15")):
        if line[0].isdigit():
            outfile.write(f"{n} {line}")

# names
with open(f"{DATA_FILE}_names.txt", "w", encoding="utf-8") as outfile:
    for n, line in enumerate(open(DATA_FILE, "r", encoding="ISO-8859-15")):
        if "MTekst 2" in str(line):
            outfile.write(f"{n}, {line[18:]}")


# In[11]:


data = pd.read_csv(
    f"{DATA_FILE}_data.txt",
    header=None,
    sep=" ",
    names=["y", "x", "z", "text"],
)

names = pd.read_csv(
    f"{DATA_FILE}_names.txt", header=None, names=["index", "name"]
)


# In[12]:


data.head()


# In[13]:


data.shape


# In[14]:


names.head()


# In[15]:


names.shape


# In[16]:


names.set_index("index", inplace=True)


# ## Merge names with coordinates

# In[17]:


names = names.reindex(range(max(data.index) + 1))


# In[18]:


# handle missing data with backfill / forward fill
names = names.bfill()
names = names.ffill()


# In[19]:


names.head()


# In[20]:


names.shape


# In[21]:


# merge names and data
data = pd.merge(data, names, left_index=True, right_index=True)


# In[22]:


data.head()


# In[23]:


data.shape


# ## Convert to geodataframe

# In[24]:


# drop duplicate entries using coordinates
data = data.drop_duplicates(["y", "x"])


# In[25]:


data.shape


# In[26]:


# convert coords from minutes to degrees
# https://gis.stackexchange.com/a/241922
data["x"] = data["x"] / 60
data["y"] = data["y"] / 60


# In[27]:


data.head()


# In[28]:


# convert to geodataframe
data = gpd.GeoDataFrame(
    data, geometry=gpd.points_from_xy(data.x, data.y, crs=4326)
)


# In[29]:


data.head()


# In[30]:


data.crs


# In[31]:


ax = data.to_crs(3857).plot(marker=".", markersize=1, figsize=(9, 9))
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


# ## Dissolve geometries

# In[32]:


# dissolve by name
data = data.dissolve("name").reset_index()


# In[33]:


data.head()


# In[34]:


data.shape


# ## Data for the Irish Sea

# In[35]:


# extent
mask = (-7, 53, -4.5, 54)


# In[36]:


data_ie = data.clip(mask)


# In[37]:


data_ie


# In[38]:


data_ie.shape


# ## Convert multipoint geometries to multilines

# In[39]:


# convert all multi points to multi lines
data_ie_1 = data_ie[
    data_ie["geometry"].astype(str).str.contains("MULTI")
].reset_index(drop=True)

data_ie_1["geometry"] = gpd.GeoSeries.from_wkt(
    "LINESTRING ("
    + data_ie_1["geometry"].astype(str).str.split("(", expand=True)[1],
    crs=4326,
)


# In[40]:


# merge multilines with remaining geometry data
data_ie = pd.concat(
    [
        data_ie[~data_ie["geometry"].astype(str).str.contains("MULTI")],
        data_ie_1,
    ]
).reset_index(drop=True)


# In[41]:


data_ie


# In[42]:


# get bounds of Irish Sea data
xmin, ymin, xmax, ymax = data_ie.to_crs(rd.CRS).total_bounds


# In[43]:


# view plotter points in the Irish Sea
ax = data.to_crs(rd.CRS).plot(alpha=0.5, figsize=(9, 9))
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


# In[44]:


# remove incorrect data lines that are obvious - Hibernia Atlantic
data_ie = data_ie.drop([6]).reset_index(drop=True)


# In[45]:


ax = (
    gpd.GeoDataFrame(geometry=data_ie.to_crs(rd.CRS).buffer(750))
    .dissolve()
    .plot(figsize=(12, 12), alpha=0.25, color="slategrey")
)
data_ie.to_crs(rd.CRS).plot(
    column="name",
    legend=True,
    ax=ax,
    cmap="jet",
    legend_kwds={"loc": "upper right"},
)
data.to_crs(rd.CRS).plot(
    ax=ax, marker="x", color="black", markersize=20, alpha=0.5
)
plt.xlim(xmin - 10000, xmax + 50000)
plt.ylim(ymin, ymax)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


# In[46]:


data_ie.to_file(os.path.join(DATA_DIR, "KIS-ORCA.gpkg"))
