#!/usr/bin/env python
# coding: utf-8

# # Kish Basin Salt Caverns

import glob
import os

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
import xarray as xr
from geocube.api.core import make_geocube
from matplotlib_scalebar.scalebar import ScaleBar

# data directory
DATA_DIR = os.path.join("data", "kish-basin")

crs = 23029

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# ## Read data layers


def read_dat_file(dat_path: str, dat_crs):
    """
    Read XYZ data layers into an Xarray dataset
    """

    gdf = {}
    # for dat_file in glob.glob(os.path.join(dat_path, "*.dat")):
    for dat_file in [
        x
        for x in glob.glob(os.path.join(dat_path, "*.dat"))
        if not "Zone" in x
    ]:
        # read each layer as individual dataframes into a dictionary
        gdf[os.path.split(dat_file)[-1][:-4]] = pd.read_fwf(
            dat_file, header=None, names=["X", "Y", "Z"]
        )

        # assign layer name to a column
        gdf[os.path.split(dat_file)[-1][:-4]]["data"] = os.path.split(
            dat_file
        )[-1][:-4]

    # find data resolution
    gdf_xr = (
        gdf[os.path.split(dat_file)[-1][:-4]].set_index(["X", "Y"]).to_xarray()
    )
    resx = gdf_xr["X"][1] - gdf_xr["X"][0]
    resy = gdf_xr["Y"][1] - gdf_xr["Y"][0]

    # combine dataframes
    gdf = pd.concat(gdf.values())

    # convert dataframe to geodataframe
    gdf["geometry"] = gpd.points_from_xy(gdf.X, gdf.Y, crs=dat_crs)
    gdf = gpd.GeoDataFrame(gdf)
    gdf.drop(columns=["X", "Y"], inplace=True)

    # convert to Xarray dataset
    ds = make_geocube(
        vector_data=gdf,
        resolution=(-abs(resy), abs(resx)),
        align=(abs(resy / 2), abs(resx / 2)),
        group_by="data",
    )

    # split variables and halite members
    ds_ = {}
    for dat in ds["data"].values:
        halite_member = dat.split(" ")[0]
        if halite_member == "Presall":
            halite_member = "Preesall"
        unit = dat.split(" ")[-1]
        zvar = dat.split("Halite ")[-1].split(" XYZ")[0]
        ds_[dat] = (
            ds.sel(data=dat)
            .assign_coords(halite=halite_member)
            .expand_dims(dim="halite")
            .drop_vars("data")
        )
        ds_[dat] = ds_[dat].rename({"Z": zvar.replace(" ", "")})
        ds_[dat][zvar.replace(" ", "")] = ds_[dat][
            zvar.replace(" ", "")
        ].assign_attrs(units=unit, long_name=zvar)

    ds = xr.combine_by_coords(ds_.values(), combine_attrs="override")

    # # keep only points corresponding to zones of interest in the dataframe
    # zones = gdf.loc[gdf["data"].str.contains("Zone")]

    # # create zones of interest polygon
    # zones = gpd.GeoDataFrame(geometry=zones.buffer(100).envelope).dissolve()

    # create extent polygon
    extent = pd.read_csv(
        os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv"), skiprows=2
    )
    extent = gpd.GeoSeries(
        shapely.geometry.Polygon(
            [
                (extent[" X"][0], extent[" Y"][0]),
                (extent[" X"][1], extent[" Y"][1]),
                (extent[" X"][2], extent[" Y"][2]),
                (extent[" X"][3], extent[" Y"][3]),
            ]
        ),
        crs=crs,
    )

    return (
        ds,
        extent,
    )  # zones


ds, extent = read_dat_file(DATA_DIR, dat_crs=crs)

xmin, ymin, xmax, ymax = extent.total_bounds


def plot_facet_maps(ds, extent, crs):
    """
    Helper function to plot facet maps of the halite layers
    """

    xmin, ymin, xmax, ymax = extent.total_bounds

    for v in ds.data_vars:
        fig = ds[v].plot.contourf(
            col="halite",
            robust=True,
            levels=15,
            cmap="jet",
            col_wrap=2,
            subplot_kws={"projection": ccrs.epsg(crs)},
            xlim=(xmin, xmax),
            ylim=(ymin, ymax),
        )
        # add a basemap
        basemap = cx.providers.CartoDB.PositronNoLabels
        for n, axis in enumerate(fig.axs.flat):
            cx.add_basemap(axis, crs=crs, source=basemap, attribution=False)
            # add attribution for basemap tiles
            if n == 2:
                axis.text(
                    xmin, ymin - 2500, basemap["attribution"], fontsize=8
                )
        fig.set_titles("{value}", weight="semibold")
        plt.show()


ds

ds.rio.crs

ds.rio.resolution()

ds.rio.bounds()

plot_facet_maps(ds, extent, crs)

# compare depths
(ds["BaseDepth"] - ds["TopDepth"]).plot.contourf(
    col="halite",
    col_wrap=2,
    levels=15,
    extend="both",
    subplot_kws={"projection": ccrs.epsg(crs)},
)
plt.show()

min(set((ds["BaseDepth"] - ds["TopDepth"]).values.flatten()))

max(set((ds["BaseDepth"] - ds["TopDepth"]).values.flatten()))

# ## Zones of interest


def zones_of_interest(ds, extent, crs, min_thickness, min_depth, max_depth):
    """
    Generate a (multi)polygon of the zones of interest by applying thickness
    and depth constraints.

    Parameters
    ----------
    ds : Xarray dataset of the halite data
    extent : extent of the data
    crs : CRS
    min_thickness : minimum halite thickness [m]
    min_depth : minimum halite top depth [m]
    max_depth : maximum halite top depth [m]

    Returns
    -------
    - A (multi)polygon geodataframe of the zones of interest
    """

    xmin, ymin, xmax, ymax = extent.total_bounds

    zdf = ds.where(
        (
            (ds.Thickness >= min_thickness)
            & (ds.TopDepth >= min_depth)
            & (ds.TopDepth <= max_depth)
        ),
        drop=True,
    )

    # zones of interest polygon
    zdf = (
        zdf.max(dim="halite")["Thickness"]
        .to_dataframe()
        .dropna()
        .reset_index()
    )
    zdf = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries(gpd.points_from_xy(zdf.x, zdf.y))
        .buffer(100)
        .envelope,
        crs=crs,
    ).dissolve()

    ax = plt.axes(projection=ccrs.epsg(crs))
    zdf.boundary.plot(ax=ax, linewidth=1, color="darkslategrey")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    cx.add_basemap(ax, source=cx.providers.CartoDB.PositronNoLabels, crs=crs)
    plt.title("Zones of interest")
    plt.tight_layout()
    plt.show()

    return zdf


# recommendations from the Hystories project
# salt thickness >= 300 m, top depth >= 1000 m and <= 1500 m
zdf = zones_of_interest(ds, extent, crs, 300, 1000, 1500)

# ## Generate potential salt cavern locations


def plot_map(data, extent, crs, var, stat):
    """
    Helper function to plot halite layer and caverns within the zones of
    interest
    """

    xmin, ymin, xmax, ymax = extent.total_bounds

    plt.figure(figsize=(12, 9))
    ax = plt.axes(projection=ccrs.epsg(crs))

    cbar_label = f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"

    if stat == "max":
        plot_data = data.max(dim="halite", skipna=True)
        cbar_label = f"Maximum {cbar_label}"
    elif stat == "min":
        plot_data = data.min(dim="halite", skipna=True)
        cbar_label = f"Minimum {cbar_label}"
    elif stat == "mean":
        plot_data = data.mean(dim="halite", skipna=True)
        cbar_label = f"Mean {cbar_label}"

    plot_data[var].plot.contourf(
        cmap="jet",
        alpha=0.65,
        robust=True,
        levels=15,
        xlim=(xmin, xmax),
        ylim=(ymin, ymax),
        cbar_kwargs={"label": cbar_label},
    )
    caverns.centroid.plot(
        ax=ax, markersize=7, color="black", label="Cavern", edgecolor="none"
    )
    cx.add_basemap(ax, crs=crs, source=cx.providers.CartoDB.Voyager)
    ax.gridlines(
        draw_labels={"bottom": "x", "left": "y"},
        alpha=0.25,
        color="darkslategrey",
    )
    ax.add_artist(
        ScaleBar(
            1,
            box_alpha=0,  # font_properties={"size": "large"},
            location="lower right",
            color="darkslategrey",
        )
    )
    plt.legend(loc="lower right", bbox_to_anchor=(1, 0.05), markerscale=1.75)
    # plt.title("Potential Kish Basin Caverns within the Zones of Interest")
    plt.title(None)
    plt.tight_layout()
    plt.show()


# ### Cavern calculations based on Caglayan *et al.* (2020)
#
# <https://doi.org/10.1016/j.ijhydene.2019.12.161>

# #### Cavern distribution function


def generate_caverns_caglayan_etal(extent, crs, diameter, separation):
    """
    Generate salt caverns using a regular square grid within the zones of
    interest based on the methodology by Caglayan et al. (2020):
    https://doi.org/10.1016/j.ijhydene.2019.12.161.
    Gridding method based on
    https://james-brennan.github.io/posts/fast_gridding_geopandas/.

    Parameters
    ----------
    extent : extent of the data
    crs : CRS
    diameter : diameter of the cavern [m]
    separation : cavern separation distance [m]

    Returns
    -------
    - A polygon geodataframe of potential caverns
    """

    xmin, ymin, xmax, ymax = extent.total_bounds

    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax + separation, separation):
        for y0 in np.arange(ymin, ymax + separation, separation):
            # bounds
            x1 = x0 - separation
            y1 = y0 + separation
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
    grid_cells = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=crs)

    # verify separation distance
    if x0 - x1 == y1 - y0:
        # generate caverns within the zones of interest
        caverns = gpd.sjoin(
            gpd.GeoDataFrame(
                geometry=grid_cells.centroid.buffer(diameter / 2)
            ).drop_duplicates(),
            zdf,
            predicate="within",
        )

        # estimates
        print("Number of potential caverns:", len(caverns))
        print(
            "Total volume:",
            "{:.2E}".format(len(caverns) * 5e5),
            f"m\N{SUPERSCRIPT THREE}",
        )
        print(
            "Estimated storage capacity:",
            "{:.2f}".format(len(caverns) * 146.418 / 1e3),
            "TWh",
        )

    else:
        print("x and y separation distances do not match!")

    return caverns


# #### Cavern capacity calculations


def cavern_capacity_caglayan_etal(
    depth, m, z, v_cavern, lhv_gas, cavern_height=120
):
    """
    t_avg : average gas temperature [K]
    depth : depth of the bottom tip of the salt cavern [m]
    cavern_height : height of the cavern [120 m]
    p_overburden : overburden / lithostatic pressure
    rho_rock : rock density [kg m-3]
    g : gravitational acceleration [9.81 m s-2]
    rho_h2 : density of hydrogen gas [kg m-3]
    z : compressibility factor
    p : pressure [Pa]
    m : molar mass of the species [kg mol-1]
    r : universal gas constant [8.314 J K-1 mol-1]
    t : temperature [K]
    m_working_gas : mass of the working gas [kg]
    v_cavern : cavern volume [m3]
    theta_safety : safety factor [70%]
    cavern_capacity : energy storage capacity of cavern [GWhH2]
    lhv_gas : lower heating value of the gas [GWhH2 kg-1]
    """

    # eq. (1)
    t_avg = 288 + 0.025 * (depth - cavern_height / 2)

    # eq. (2)
    # p_overburden = rho_rock * g * (depth - cavern_height)
    p_overburden = rho_rock * 9.81 * (depth - cavern_height)

    # eq. (3)
    # rho_h2 = (p * m) / (z * r * t)
    rho_h2_max = (p_overburden * 0.8 * m) / (z * 8.314 * t_avg)
    rho_h2_min = (p_overburden * 0.24 * m) / (z * 8.314 * t_avg)

    # eq. (4)
    # m_working_gas = (rho_h2_max - rho_h2_min) * v_cavern * theta_safety
    m_working_gas = (rho_h2_max - rho_h2_min) * v_cavern * 0.7

    # eq. (5)
    cavern_capacity = m_working_gas * lhv_gas


# #### Generate caverns

# max 80 m diameter based on Hystories
# separation distance of 4 times the diameter as recommended by Caglayan et al.
caverns = generate_caverns_caglayan_etal(extent, crs, 80, 80 * 4)

# 85 m diameter and 330 m separation distance as used by HYSS
caverns = generate_caverns_caglayan_etal(extent, crs, 85, 330)

plot_map(ds, extent, crs, "Thickness", "max")

# ### Cavern calculations based on Williams *et al.* (2022)
#
# <https://doi.org/10.1016/j.est.2022.105109>

# #### Cavern distribution function


def generate_caverns_williams_etal(extent, crs, diameter, separation):
    """
    Generate caverns in a regular hexagonal grid as proposed by Williams
    et al. (2022): https://doi.org/10.1016/j.est.2022.105109.
    Hexagonal gridding method based on
    https://sabrinadchan.github.io/data-blog/building-a-hexagonal-cartogram.html.

    Parameters
    ----------
    extent : extent of the data
    crs : CRS
    diameter : diameter of the cavern [m]
    separation : cavern separation distance [m]

    Returns
    -------
    - A polygon geodataframe of potential caverns
    """

    xmin, ymin, xmax, ymax = extent.total_bounds

    a = np.sin(np.pi / 3)
    cols = np.arange(np.floor(xmin), np.ceil(xmax), 3 * separation)
    rows = np.arange(np.floor(ymin) / a, np.ceil(ymax) / a, separation)

    hexagons = []
    caverns = []
    for x in cols:
        for i, y in enumerate(rows):
            if i % 2 == 0:
                x0 = x
            else:
                x0 = x + 1.5 * separation

            hexagons.append(
                shapely.geometry.Polygon(
                    [
                        (x0, y * a),
                        (x0 + separation, y * a),
                        (x0 + (1.5 * separation), (y + separation) * a),
                        (x0 + separation, (y + (2 * separation)) * a),
                        (x0, (y + (2 * separation)) * a),
                        (x0 - (0.5 * separation), (y + separation) * a),
                    ]
                )
            )

            caverns.append(shapely.geometry.Point(x0, y * a))
            caverns.append(shapely.geometry.Point(x0 + separation, y * a))
            caverns.append(
                shapely.geometry.Point(
                    x0 + (1.5 * separation), (y + separation) * a
                )
            )
            caverns.append(
                shapely.geometry.Point(
                    x0 + separation, (y + (2 * separation)) * a
                )
            )
            caverns.append(
                shapely.geometry.Point(x0, (y + (2 * separation)) * a)
            )
            caverns.append(
                shapely.geometry.Point(
                    x0 - (0.5 * separation), (y + separation) * a
                )
            )

    hexagons = gpd.GeoDataFrame(geometry=hexagons, crs=crs)

    # generate caverns using hexagon vertices and centroids
    caverns = pd.concat(
        [
            gpd.GeoDataFrame(geometry=caverns, crs=crs),
            gpd.GeoDataFrame(geometry=hexagons.centroid, crs=crs),
        ]
    ).drop_duplicates()
    caverns = gpd.GeoDataFrame(geometry=caverns.buffer(diameter / 2))

    # clip caverns to the zones of interest
    caverns = gpd.sjoin(caverns, zdf, predicate="within")

    # estimates
    print("Number of potential caverns:", len(caverns))
    print(
        "Total volume:",
        "{:.2E}".format(len(caverns) * 5e5),
        f"m\N{SUPERSCRIPT THREE}",
    )
    print(
        "Estimated storage capacity:",
        "{:.2f}".format(len(caverns) * 105.074 / 1e3),
        "TWh",
    )

    return caverns


# max 80 m diameter based on Hystories
# separation distance of 4 times the diameter as recommended by Caglayan et al.
caverns = generate_caverns_williams_etal(extent, crs, 80, 80 * 4)

# 85 m diameter and 330 m separation distance as used by HYSS
caverns = generate_caverns_williams_etal(extent, crs, 85, 330)

plot_map(ds, extent, crs, "Thickness", "max")
