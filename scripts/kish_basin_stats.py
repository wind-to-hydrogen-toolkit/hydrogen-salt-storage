#!/usr/bin/env python
# coding: utf-8

# # Kish Basin statistics

import os

import cartopy.crs as ccrs
import contextily as cx
import matplotlib.pyplot as plt
import seaborn as sns
from cartopy.mpl.ticker import LongitudeFormatter
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import data as rd

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# ## Read data layers

ds, extent = rd.kish_basin_data_depth_adjusted(
    dat_path=os.path.join("data", "kish-basin"),
    bathymetry_path=os.path.join("data", "bathymetry"),
)

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


plot_facet_maps(ds, extent, rd.CRS)

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
        vars=["Thickness", "TopDepthSeabed", "BaseDepthSeabed"],
    )
    plt.show()

    # histograms
    _, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
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
        x="TopDepthSeabed",
        hue="halite",
        ax=axes[1],
        palette="rocket_r",
        multiple="fill",
        bins=100,
    )
    plt.tight_layout()
    plt.show()

    # box plots
    _, axes = plt.subplots(1, 3, figsize=(8, 4.5))
    sns.boxplot(
        dat_df.reset_index(),
        y="Thickness",
        hue="halite",
        palette="flare",
        ax=axes[0],
    )
    sns.boxplot(
        dat_df.reset_index(),
        y="TopDepthSeabed",
        hue="halite",
        palette="flare",
        ax=axes[1],
        legend=False,
    )
    sns.boxplot(
        dat_df.reset_index(),
        y="BaseDepthSeabed",
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
shape = rd.halite_shape(dat_xr=ds)

f"Surface area: {shape.area[0]:.2E} m\N{SUPERSCRIPT TWO}"


def plot_facet_maps_distr(
    dat_xr,
    dat_extent,
    dat_crs,
    v,
    levels,
    label,
    scalebar=True,
    attribution=True,
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
            "pad": 0.09,
            "extendfrac": 0.2,
            "label": label,
            "format": lambda x, _: f"{x:,.0f}",
        },
        col_wrap=2,
    )

    # add a basemap
    basemap = cx.providers.CartoDB.PositronNoLabels
    for n, axis in enumerate(f.axs.flat):
        cx.add_basemap(axis, crs=dat_crs, source=basemap, attribution=False)
        # tick labels and attribution for basemap tiles
        if n in (0, 2):
            axis.gridlines(draw_labels={"left": "y"}, color="none")
        if n in (2, 3):
            axis.gridlines(
                draw_labels={"bottom": "x"},
                color="none",
                rotate_labels=90,
                xpadding=20,
                xformatter=LongitudeFormatter(number_format=".2f"),
            )
        if n == 3:
            if scalebar:
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
            if attribution:
                axis.text(
                    xmin_ + 1000,
                    ymin_ + 1000,
                    basemap["attribution"],
                    fontsize=9,
                )
    f.set_titles("{value}", weight="semibold")
    plt.show()


plot_facet_maps_distr(
    ds,
    extent,
    rd.CRS,
    "TopDepthSeabed",
    [500 - 80, 1000 - 80, 1500 - 80, 2000 - 80],
    "Top Depth [m]",
    scalebar=False,
)

plot_facet_maps_distr(
    ds,
    extent,
    rd.CRS,
    "Thickness",
    [85 + 90, 155 + 90, 311 + 90],
    "Thickness [m]",
    attribution=False,
)

# compare depths
(ds["BaseDepth"] - ds["TopDepth"]).plot(
    col="halite",
    col_wrap=2,
    extend="both",
    subplot_kws={"projection": ccrs.epsg(rd.CRS)},
)
plt.show()

min(set((ds["BaseDepth"] - ds["TopDepth"]).values.flatten()))

max(set((ds["BaseDepth"] - ds["TopDepth"]).values.flatten()))
