#!/usr/bin/env python
# coding: utf-8

# # Ireland boundary
#
# Provinces - OSi National Statutory Boundaries - 2019 - Ungeneralised
#
# - <https://data.gov.ie/dataset/provinces-national-statutory-boundaries-2019>
# - <https://data-osi.opendata.arcgis.com/maps/osi::provinces-national-statutory-boundaries-2019>

# In[1]:


import os
from zipfile import ZipFile

import matplotlib.pyplot as plt

from h2ss import data as rd

# In[2]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0


# In[3]:


URL = (
    "https://data-osi.opendata.arcgis.com/datasets/"
    "559bc3300384413aa0fe93f0772cb7f1_0.zip?"
    "outSR=%7B%22latestWkid%22%3A2157%2C%22wkid%22%3A2157%7D"
)
DATA_DIR = os.path.join("data", "boundaries")
FILE_NAME = "osi-provinces-ungeneralised-2019.zip"
DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)


# In[4]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# In[5]:


ZipFile(DATA_FILE).namelist()


# In[6]:


data = rd.read_shapefile_from_zip(data_path=DATA_FILE)


# In[7]:


data


# In[8]:


data.crs


# In[9]:


data.shape


# In[10]:


data = data.dissolve()


# In[11]:


data.bounds


# In[12]:


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
