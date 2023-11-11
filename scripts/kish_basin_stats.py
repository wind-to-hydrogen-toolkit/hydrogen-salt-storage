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
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
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

# surface area
shape = fns.halite_shape(ds, CRS)

f"Surface area: {shape.area[0]:.2E} m\N{SUPERSCRIPT TWO}"


def plot_facet_maps_distr(dat_xr, dat_extent, dat_crs, v, levels):
    """
    Helper function to plot facet maps of the halite layers

    Parameters
    ----------
    dat_xr : Xarray dataset of the halite data
    dat_extent : extent of the data
    dat_crs : EPSG CRS
    """

    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds

    f = dat_xr[v].plot.contourf(
        col="halite",
        robust=True,
        levels=levels,
        cmap=sns.color_palette("flare", as_cmap=True),
        # col_wrap=2,
        figsize=(12, 6),
        subplot_kws={"projection": ccrs.epsg(dat_crs)},
        xlim=(xmin_, xmax_),
        ylim=(ymin_, ymax_),
        cbar_kwargs={
            "location": "bottom",
            "aspect": 35,
            "shrink": 0.65,
            "pad": 0.05,
            "extendfrac": 0.2,
            "label": dat_xr[v].attrs["long_name"] + " [m]",
        },
    )
    # add a basemap
    basemap = cx.providers.CartoDB.PositronNoLabels
    for n, axis in enumerate(f.axs.flat):
        cx.add_basemap(axis, crs=dat_crs, source=basemap, attribution=False)
        # add attribution for basemap tiles
        if n == 0:
            axis.text(xmin_, ymin_ - 2500, basemap["attribution"], fontsize=8)
    f.set_titles("{value}", weight="semibold")
    plt.show()


plot_facet_maps_distr(
    ds, extent, CRS, "TopDepth", [500 - 80, 1000 - 80, 1500 - 80, 2000 - 80]
)

plot_facet_maps_distr(
    ds, extent, CRS, "Thickness", [85 + 90, 155 + 90, 311 + 90]
)

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


def plot_zones_map(zdf, dat_extent, dat_crs):
    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds

    ax = plt.axes(projection=ccrs.epsg(dat_crs))
    zdf.boundary.plot(ax=ax, linewidth=1, color="darkslategrey")
    plt.xlim(xmin_, xmax_)
    plt.ylim(ymin_, ymax_)
    cx.add_basemap(
        ax, source=cx.providers.CartoDB.PositronNoLabels, crs=dat_crs
    )
    plt.title("Zones of interest")
    plt.tight_layout()
    plt.show()


# ### Case 1
#
# height = 311 m, 1,000 m <= depth <= 1,500 m, diameter = 80 m,
# separation = 320 m

zones_1, zds_1 = fns.zones_of_interest(
    ds, CRS, {"height": 311, "min_depth": 1000, "max_depth": 1500}
)

plot_zones_map(zones_1, extent, CRS)

# ### Case 2
#
# height = 155 m, 1000 m <= depth <= 1,500 m, diameter = 80 m,
# separation = 320 m

zones_2, zds_2 = fns.zones_of_interest(
    ds, CRS, {"height": 155, "min_depth": 750, "max_depth": 1750}
)

plot_zones_map(zones_2, extent, CRS)

# ### Case 3
#
# height = 85 m, 500 m <= depth <= 2,000 m, diameter = 80 m,
# separation = 320 m

zones_3, zds_3 = fns.zones_of_interest(
    ds, CRS, {"height": 85, "min_depth": 500, "max_depth": 2000}
)

plot_zones_map(zones_3, extent, CRS)

# ## Zones of interest stats

# ### Case 1

zdf_1 = make_stats_plots(zds_1)

zdf_1.describe()

# ### Case 2

zdf_2 = make_stats_plots(zds_2)

zdf_2.describe()

# ### Case 3

zdf_3 = make_stats_plots(zds_3)

zdf_3.describe()

# ## Sensitivity analysis
#
# Sensitivity of maximum halite area (i.e. zones of interest) available for
# cavern construction (without constraints) to depth and thickness


def sensitivity(zones_gdf, height_base, min_depth_base, max_depth_base):
    sdf = {}

    # sensitivity to height
    heights = [85 + 226 / 50 * n for n in range(51)]
    area_max = []

    for t in heights:
        area_max.append(
            fns.zones_of_interest(
                ds,
                CRS,
                {
                    "height": t,
                    "min_depth": min_depth_base,
                    "max_depth": max_depth_base,
                },
            )[0].area[0]
        )

    sdf["height"] = pd.DataFrame({"max_area": area_max, "height": heights})

    # percentage change
    sdf["height"]["diff_area"] = (
        (sdf["height"]["max_area"] - zones_gdf.area[0])
        / zones_gdf.area[0]
        * 100
    )

    # sensitivity to minimum depth
    depth_min = [500 + 10 * n for n in range(61)]
    area_max = []

    for t in depth_min:
        area_max.append(
            fns.zones_of_interest(
                ds,
                CRS,
                {
                    "height": height_base,
                    "min_depth": t,
                    "max_depth": max_depth_base,
                },
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

    # sensitivity to maximum depth
    depth_max = [1400 + 10 * n for n in range(61)]
    area_max = []

    for t in depth_max:
        area_max.append(
            fns.zones_of_interest(
                ds,
                CRS,
                {
                    "height": height_base,
                    "min_depth": min_depth_base,
                    "max_depth": t,
                },
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


sdf = sensitivity(zones_2, 155, 750, 1750)

for key in sdf.keys():
    print(key)
    print(sdf[key].describe())

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for (n, key), c in zip(
    enumerate(sdf.keys()), ["royalblue", "crimson", "seagreen"]
):
    sns.regplot(
        sdf[key],
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
    enumerate(sdf.keys()), ["royalblue", "crimson", "seagreen"]
):
    sns.regplot(
        sdf[key],
        x=key,
        y="diff_area",
        logx=True,
        ax=axes[n],
        marker=".",
        color=c,
    )
plt.ylim(-110, 110)
plt.tight_layout()
plt.show()
