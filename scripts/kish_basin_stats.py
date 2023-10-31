#!/usr/bin/env python
# coding: utf-8

# # Kish Basin Statistics

import os

import cartopy.crs as ccrs
import contextily as cx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src import functions as fns

# data directory
DATA_DIR = os.path.join("data", "kish-basin")

CRS = 23029

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# ## Read data layers

ds, extent = fns.read_dat_file(DATA_DIR, CRS)

ds

ds.rio.crs

ds.rio.resolution()

ds.rio.bounds()


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


plot_facet_maps(ds, extent, CRS)

# ## Stats

df = ds.to_dataframe()[list(ds.data_vars)]

df.describe()

sns.pairplot(
    df.reset_index(),
    palette="flare",
    hue="halite",
    plot_kws={"alpha": 0.5},
    vars=["Thickness", "TopDepth", "BaseDepth"],
)
plt.show()

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sns.histplot(
    df.reset_index(),
    x="Thickness",
    hue="halite",
    ax=axes[0],
    palette="rocket_r",
    multiple="fill",
    bins=100,
    legend=False,
)
sns.histplot(
    df.reset_index(),
    x="TopDepth",
    hue="halite",
    ax=axes[1],
    palette="rocket_r",
    multiple="fill",
    bins=100,
)
plt.tight_layout()
plt.show()

plt.subplots(figsize=(12, 5))
sns.boxplot(
    pd.melt(df.reset_index().drop(columns=["x", "y"]), id_vars=["halite"]),
    x="variable",
    y="value",
    hue="halite",
    palette="flare",
)
plt.show()

# compare depths
(ds["BaseDepth"] - ds["TopDepth"]).plot(
    col="halite",
    col_wrap=2,
    extend="both",
    subplot_kws={"projection": ccrs.epsg(CRS)},
)
plt.show()

min(set((ds["BaseDepth"] - ds["TopDepth"]).values.flatten()))

max(set((ds["BaseDepth"] - ds["TopDepth"]).values.flatten()))

# ## Zones of interest

# numbers used in HYSS calculations
# thickness >= 300 m, 1000 m <= depth <= 1500 m, diameter = 85 m
# separation = 330 m
zones, zones_ds = fns.zones_of_interest(
    ds,
    extent,
    CRS,
    {"min_thickness": 300, "min_depth": 1000, "max_depth": 1500},
)

# ## Zones of interest stats

zones_df = zones_ds.to_dataframe()[list(ds.data_vars)]

zones_df.describe()

sns.pairplot(
    zones_df.reset_index(),
    palette="flare",
    hue="halite",
    plot_kws={"alpha": 0.5},
    vars=["Thickness", "TopDepth", "BaseDepth"],
)
plt.show()

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sns.histplot(
    zones_df.reset_index(),
    x="Thickness",
    hue="halite",
    ax=axes[0],
    palette="rocket_r",
    multiple="fill",
    bins=100,
    legend=False,
)
sns.histplot(
    zones_df.reset_index(),
    x="TopDepth",
    hue="halite",
    ax=axes[1],
    palette="rocket_r",
    multiple="fill",
    bins=100,
)
plt.tight_layout()
plt.show()

plt.subplots(figsize=(12, 5))
sns.boxplot(
    pd.melt(
        zones_df.reset_index().drop(columns=["x", "y"]), id_vars=["halite"]
    ),
    x="variable",
    y="value",
    hue="halite",
    palette="flare",
)
plt.show()

# ## Sensitivity analysis
#
# Sensitivity of maximum halite area (i.e. zones of interest) available for
# cavern construction (without constraints) to depth and thickness

# sensitivity to minimum thickness
thickness_min = [200 + 5 * n for n in range(41)]
area_max = []

for t in thickness_min:
    area_max.append(
        fns.zones_of_interest(
            ds,
            extent,
            CRS,
            {"min_thickness": t, "min_depth": 1000, "max_depth": 1500},
            display_map=False,
        )[0].area[0]
    )

sdf = pd.DataFrame({"max_area": area_max, "min_thickness": thickness_min})

# percentage change
sdf["diff_area"] = (sdf["max_area"] - zones.area[0]) / zones.area[0] * 100

sdf.describe()

# sensitivity to minimum depth
depth_min = [500 + 12.5 * n for n in range(49)]
area_max = []

for t in depth_min:
    area_max.append(
        fns.zones_of_interest(
            ds,
            extent,
            CRS,
            {"min_thickness": 300, "min_depth": t, "max_depth": 1500},
            display_map=False,
        )[0].area[0]
    )

sdf1 = pd.DataFrame({"max_area": area_max, "min_depth": depth_min})

# percentage change
sdf1["diff_area"] = (sdf1["max_area"] - zones.area[0]) / zones.area[0] * 100

sdf1.describe()

# sensitivity to maximum depth
depth_max = [1400 + 12.5 * n for n in range(57)]
area_max = []

for t in depth_max:
    area_max.append(
        fns.zones_of_interest(
            ds,
            extent,
            CRS,
            {"min_thickness": 300, "min_depth": 1000, "max_depth": t},
            display_map=False,
        )[0].area[0]
    )

sdf2 = pd.DataFrame({"max_area": area_max, "max_depth": depth_max})

# percentage change
sdf2["diff_area"] = (sdf2["max_area"] - zones.area[0]) / zones.area[0] * 100

sdf2.describe()

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
sns.regplot(
    sdf, x="min_thickness", y="max_area", logx=True, ax=axes[0], marker="."
)
sns.regplot(
    sdf1,
    x="min_depth",
    y="max_area",
    logx=True,
    ax=axes[1],
    color="crimson",
    marker=".",
)
sns.regplot(
    sdf2,
    x="max_depth",
    y="max_area",
    logx=True,
    ax=axes[2],
    color="seagreen",
    marker=".",
)
plt.tight_layout()
plt.show()

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
sns.regplot(
    sdf, x="min_thickness", y="diff_area", logx=True, ax=axes[0], marker="."
)
sns.regplot(
    sdf1,
    x="min_depth",
    y="diff_area",
    logx=True,
    ax=axes[1],
    color="crimson",
    marker=".",
)
sns.regplot(
    sdf2,
    x="max_depth",
    y="diff_area",
    logx=True,
    ax=axes[2],
    color="seagreen",
    marker=".",
)
plt.tight_layout()
plt.show()
