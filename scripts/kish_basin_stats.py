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


def make_stats_plots(dat_xr):
    """
    Statistical plots for the halite Xarray dataset

    Parameters
    ----------
    dat_xr : Xarray dataset of the halite data

    Returns
    -------
    - The dataset converted into a dataframe
    """

    # convert to dataframe
    dat_df = dat_xr.to_dataframe()[list(dat_xr.data_vars)]

    # pairwise comparison of variables
    sns.pairplot(
        dat_df.reset_index(),
        palette="flare",
        hue="halite",
        plot_kws={"alpha": 0.5},
        vars=["Thickness", "TopDepth", "BaseDepth"],
    )
    plt.show()

    # histograms
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.histplot(
        dat_df.reset_index(),
        x="Thickness",
        hue="halite",
        ax=axes[0],
        palette="rocket_r",
        multiple="fill",
        bins=100,
        legend=False,
    )
    sns.histplot(
        dat_df.reset_index(),
        x="TopDepth",
        hue="halite",
        ax=axes[1],
        palette="rocket_r",
        multiple="fill",
        bins=100,
    )
    plt.tight_layout()
    plt.show()

    # box plots
    fig, axes = plt.subplots(1, 3, figsize=(8, 4.5))
    sns.boxplot(
        dat_df.reset_index(),
        y="Thickness",
        hue="halite",
        palette="flare",
        ax=axes[0],
    )
    sns.boxplot(
        dat_df.reset_index(),
        y="TopDepth",
        hue="halite",
        palette="flare",
        ax=axes[1],
        legend=False,
    )
    sns.boxplot(
        dat_df.reset_index(),
        y="BaseDepth",
        hue="halite",
        palette="flare",
        ax=axes[2],
        legend=False,
    )
    plt.tight_layout()
    plt.show()

    return dat_df


df = make_stats_plots(ds)

df.describe()

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

# ### With max depth

# numbers used in HYSS calculations
# thickness >= 300 m, 1000 m <= depth <= 1500 m, diameter = 85 m
# separation = 330 m
zones_max_depth, zds_max_depth = fns.zones_of_interest(
    ds,
    extent,
    CRS,
    {"min_thickness": 300, "min_depth": 1000, "max_depth": 1500},
)

# ### Without max depth

# numbers used in HYSS calculations
# thickness >= 300 m, 1000 m <= depth <= 1500 m, diameter = 85 m
# separation = 330 m
zones_min_depth, zds_min_depth = fns.zones_of_interest(
    ds, extent, CRS, {"min_thickness": 300, "min_depth": 1000}
)

# ## Zones of interest stats

# ### With max depth

zdf_max_depth = make_stats_plots(zds_max_depth)

zdf_max_depth.describe()

# ### Without max depth

zdf_min_depth = make_stats_plots(zds_min_depth)

zdf_min_depth.describe()

# ## Sensitivity analysis
#
# Sensitivity of maximum halite area (i.e. zones of interest) available for
# cavern construction (without constraints) to depth and thickness


def sensitivity(zones_gdf, include_max_depth=True):
    sdf = {}

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

    sdf["min_thickness"] = pd.DataFrame(
        {"max_area": area_max, "min_thickness": thickness_min}
    )

    # percentage change
    sdf["min_thickness"]["diff_area"] = (
        (sdf["min_thickness"]["max_area"] - zones_gdf.area[0])
        / zones_gdf.area[0]
        * 100
    )

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

    sdf["min_depth"] = pd.DataFrame(
        {"max_area": area_max, "min_depth": depth_min}
    )

    # percentage change
    sdf["min_depth"]["diff_area"] = (
        (sdf["min_depth"]["max_area"] - zones_gdf.area[0])
        / zones_gdf.area[0]
        * 100
    )

    if include_max_depth:
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

        sdf["max_depth"] = pd.DataFrame(
            {"max_area": area_max, "max_depth": depth_max}
        )

        # percentage change
        sdf["max_depth"]["diff_area"] = (
            (sdf["max_depth"]["max_area"] - zones_gdf.area[0])
            / zones_gdf.area[0]
            * 100
        )

    return sdf


# ### With max depth

sdf_max_depth = sensitivity(zones_max_depth)

for key in sdf_max_depth.keys():
    print(key)
    print(sdf_max_depth[key].describe())

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for (n, key), c in zip(
    enumerate(sdf_max_depth.keys()), ["royalblue", "crimson", "seagreen"]
):
    sns.regplot(
        sdf_max_depth[key],
        x=key,
        y="max_area",
        logx=True,
        ax=axes[n],
        marker=".",
        color=c,
    )
plt.tight_layout()
plt.show()

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
for (n, key), c in zip(
    enumerate(sdf_max_depth.keys()), ["royalblue", "crimson", "seagreen"]
):
    sns.regplot(
        sdf_max_depth[key],
        x=key,
        y="diff_area",
        logx=True,
        ax=axes[n],
        marker=".",
        color=c,
    )
plt.tight_layout()
plt.show()

# ### Without max depth

sdf_min_depth = sensitivity(zones_min_depth, include_max_depth=False)

for key in sdf_min_depth.keys():
    print(key)
    print(sdf_min_depth[key].describe())

fig, axes = plt.subplots(1, 2, figsize=(8, 4))
for (n, key), c in zip(
    enumerate(sdf_min_depth.keys()), ["royalblue", "crimson", "seagreen"]
):
    sns.regplot(
        sdf_min_depth[key],
        x=key,
        y="max_area",
        logx=True,
        ax=axes[n],
        marker=".",
        color=c,
    )
plt.tight_layout()
plt.show()

fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
for (n, key), c in zip(
    enumerate(sdf_min_depth.keys()), ["royalblue", "crimson", "seagreen"]
):
    sns.regplot(
        sdf_min_depth[key],
        x=key,
        y="diff_area",
        logx=True,
        ax=axes[n],
        marker=".",
        color=c,
    )
plt.tight_layout()
plt.show()
