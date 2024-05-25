#!/usr/bin/env python
# coding: utf-8

# # Kish Basin statistics

# In[1]:


import os

import cartopy.crs as ccrs
import contextily as cx
import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.patches as mpatches
import seaborn as sns
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import data as rd
from h2ss import functions as fns

# In[2]:


# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[3]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0
plt.rcParams["xtick.minor.size"] = 0
plt.rcParams["ytick.minor.size"] = 0


# ## Read data layers

# In[4]:


ds, extent = rd.kish_basin_data_depth_adjusted(
    dat_path=os.path.join("data", "kish-basin"),
    bathymetry_path=os.path.join("data", "bathymetry"),
)


# In[5]:


ds = fns.net_to_gross(ds)


# In[6]:


ds


# In[7]:


ds.rio.crs


# In[8]:


ds.rio.resolution()


# In[9]:


ds.rio.bounds()


# In[25]:


def plot_facet_maps(dat_xr, dat_extent, dat_crs):
    """
    Helper function to plot facet maps of the halite layers

    Parameters
    ----------
    dat_xr : Xarray dataset of the halite data
    dat_extent : extent of the data
    dat_crs : EPSG CRS
    """

    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds

    for v in dat_xr.data_vars:
        f = dat_xr[v].plot.contourf(
            col="halite",
            robust=True,
            levels=15,
            cmap="jet",
            col_wrap=2,
            subplot_kws={"projection": ccrs.epsg(dat_crs)},
            xlim=(xmin_, xmax_),
            ylim=(ymin_, ymax_),
            cbar_kwargs={"label": v},
        )
        # add a basemap
        basemap = cx.providers.CartoDB.PositronNoLabels
        for n, axis in enumerate(f.axs.flat):
            cx.add_basemap(
                axis, crs=dat_crs, source=basemap, attribution=False
            )
            # add attribution for basemap tiles
            if n == 2:
                axis.text(
                    xmin_, ymin_ - 2500, basemap["attribution"], fontsize=8
                )
        f.set_titles("{value}", weight="semibold")
        plt.show()


# In[26]:


plot_facet_maps(ds, extent, rd.CRS)


# ## Stats

# In[12]:


df = ds.to_dataframe()[list(ds.data_vars)]


# In[13]:


df.describe()


# In[14]:


# max total thickness
ds.sum(dim="halite").to_dataframe()[["Thickness", "ThicknessNet"]].max()


# In[15]:


# surface area
shape = rd.halite_shape(dat_xr=ds)


# In[16]:


print(f"Surface area: {shape.area[0]:.4E} m\N{SUPERSCRIPT TWO}")


# In[47]:


def plot_facet_maps_distr(
    dat_xr,
    dat_extent,
    dat_crs,
    v,
    levels,
    label,
):
    """
    Helper function to plot facet maps of the halite layers

    Parameters
    ----------
    dat_xr : Xarray dataset of the halite data
    dat_extent : extent of the data
    dat_crs : EPSG CRS
    """

    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds
    if levels:
        levels = sorted(levels)
    # legend_handles = []

    f = dat_xr[v].plot.contourf(
        col="halite",
        robust=True,
        levels=levels,
        cmap=sns.color_palette("flare", as_cmap=True),
        figsize=(6, 8.5),
        subplot_kws={"projection": ccrs.epsg(dat_crs)},
        xlim=(xmin_, xmax_),
        ylim=(ymin_, ymax_),
        cbar_kwargs={
            "location": "bottom",
            "aspect": 20,
            "shrink": 0.8,
            "pad": 0.07,
            "extendfrac": 0.2,
            "label": label,
            "format": lambda x, _: f"{x:,.0f}",
        },
        # add_colorbar=False,
        col_wrap=2,
    )

    # colours = [int(n * 255 / (len(levels) - 1)) for n in range(len(levels))]
    # for n, (level, colour) in enumerate(zip(levels, colours)):
    #     if n == 0:
    #         legend_handles.append(
    #             mpatches.Patch(
    #                 facecolor=sns.color_palette("flare", 256)[colour],
    #                 label=f"< {level}"
    #             )
    #         )
    #     elif n == len(levels) - 1:
    #         legend_handles.append(
    #             mpatches.Patch(
    #                 facecolor=sns.color_palette("flare", 256)[colour],
    #                 label=f"> {levels[n - 1]}"
    #             )
    #         )
    #     else:
    #         legend_handles.append(
    #             mpatches.Patch(
    #                 facecolor=sns.color_palette("flare", 256)[colour],
    #                 label=f"{levels[n - 1]}–{level}"
    #             )
    #         )

    # add a basemap
    basemap = cx.providers.CartoDB.PositronNoLabels
    for n, axis in enumerate(f.axs.flat):
        cx.add_basemap(axis, crs=dat_crs, source=basemap, attribution=False)
        # tick labels and attribution for basemap tiles
        if n in (0, 2):
            axis.gridlines(
                draw_labels={"left": "y"},
                color="none",
                yformatter=LatitudeFormatter(auto_hide=False, dms=True),
                ylabel_style={"rotation": 89.9},
                ylocs=[53.2, 53.4],
            )
        if n in (2, 3):
            axis.gridlines(
                draw_labels={"bottom": "x"},
                color="none",
                xformatter=LongitudeFormatter(auto_hide=False, dms=True),
                xlocs=[-6.15, -5.85, -5.55],
            )
        axis.gridlines(color="lightslategrey", alpha=0.25)
        if n == 3:
            axis.add_artist(
                ScaleBar(
                    1,
                    box_alpha=0,
                    location="lower right",
                    color="darkslategrey",
                    width_fraction=0.015,
                )
            )
        if n == 2:
            axis.text(
                xmin_ + 1000,
                ymin_ + 1000,
                basemap["attribution"],
                fontsize=8,
            )
    f.set_titles("{value}", weight="semibold")
    # plt.legend(
    #     bbox_to_anchor=(0.875, -0.12),
    #     title=label,
    #     handles=legend_handles,
    #     fontsize=11.5,
    #     ncols=3
    # )
    # plt.savefig(
    #     os.path.join("graphics", f"fig_facet_{v.lower()}.jpg"),
    #     format="jpg",
    #     dpi=600,
    # )
    plt.show()


# In[48]:


plot_facet_maps_distr(
    ds,
    extent,
    rd.CRS,
    "ThicknessNet",
    [85 + 35 * n + 90 for n in range(7)],
    "Net thickness [m]",
)


# In[49]:


plot_facet_maps_distr(
    ds,
    extent,
    rd.CRS,
    "TopDepthSeabed",
    [500 - 80, 1000 - 80, 1500 - 80, 2000 - 80],
    "Top depth [m]",
)


# In[20]:


# compare depths
(ds["BaseDepth"] - ds["TopDepth"]).plot(
    col="halite",
    col_wrap=2,
    extend="both",
    subplot_kws={"projection": ccrs.epsg(rd.CRS)},
)
plt.show()


# In[21]:


val = (ds["BaseDepth"] - ds["TopDepth"]).values.flatten()
min(val[~np.isnan(val)])


# In[22]:


max(val[~np.isnan(val)])
