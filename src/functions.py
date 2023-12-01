"""
Main functions used in the Jupyter Notebooks
"""

import glob
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr
from geocube.api.core import make_geocube


def read_dat_file(
    dat_path: str, dat_crs: int = 23029
) -> tuple[xr.Dataset, gpd.GeoSeries]:
    """
    Read XYZ data layers into an Xarray dataset

    Parameters
    ----------
    dat_path : path to the .dat files
    dat_crs : EPSG CRS

    Returns
    -------
    - Xarray dataset
    - GeoPandas GeoSeries of the extent
    """

    gdf = {}
    # for dat_file in glob.glob(os.path.join(dat_path, "*.dat")):
    for dat_file in [
        x
        for x in glob.glob(os.path.join(dat_path, "*.dat"))
        if "Zone" not in x
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
    gdf_xr = gdf[list(gdf.keys())[0]].set_index(["X", "Y"]).to_xarray()
    resx = gdf_xr["X"][1] - gdf_xr["X"][0]
    resy = gdf_xr["Y"][1] - gdf_xr["Y"][0]

    # combine dataframes
    gdf = pd.concat(gdf.values())

    # convert dataframe to geodataframe
    gdf["geometry"] = gpd.points_from_xy(gdf.X, gdf.Y, crs=dat_crs)
    gdf = gpd.GeoDataFrame(gdf)
    gdf.drop(columns=["X", "Y"], inplace=True)

    # convert to Xarray dataset
    xds = make_geocube(
        vector_data=gdf,
        resolution=(-abs(resy), abs(resx)),
        align=(abs(resy / 2), abs(resx / 2)),
        group_by="data",
    )

    # split variables and halite members
    xds_ = {}
    for d in xds["data"].values:
        halite_member = d.split(" ")[0]
        # fix halite names
        if halite_member == "Presall":
            halite_member = "Preesall"
        elif halite_member == "Flyde":
            halite_member = "Fylde"
        unit = d.split(" ")[-1]
        zvar = d.split("Halite ")[-1].split(" XYZ")[0]
        xds_[d] = (
            xds.sel(data=d)
            .assign_coords(halite=halite_member)
            .expand_dims(dim="halite")
            .drop_vars("data")
        )
        xds_[d] = xds_[d].rename({"Z": zvar.replace(" ", "")})
        xds_[d][zvar.replace(" ", "")] = xds_[d][
            zvar.replace(" ", "")
        ].assign_attrs(units=unit, long_name=zvar)

    xds = xr.combine_by_coords(xds_.values(), combine_attrs="override")

    # # keep only points corresponding to zones of interest in the dataframe
    # zones = gdf.loc[gdf["data"].str.contains("Zone")]

    # # create zones of interest polygon
    # zones = gpd.GeoDataFrame(geometry=zones.buffer(100).envelope).dissolve()

    # create extent polygon
    dat_extent = pd.read_csv(
        os.path.join(dat_path, "Kish GIS Map Extent - Square.csv"), skiprows=2
    )
    dat_extent = gpd.GeoSeries(
        shapely.geometry.Polygon(
            [
                (dat_extent[" X"][0], dat_extent[" Y"][0]),
                (dat_extent[" X"][1], dat_extent[" Y"][1]),
                (dat_extent[" X"][2], dat_extent[" Y"][2]),
                (dat_extent[" X"][3], dat_extent[" Y"][3]),
            ]
        ),
        crs=dat_crs,
    )

    return xds, dat_extent  # zones


def halite_shape(
    dat_xr: xr.Dataset, dat_crs: int = 23029, halite: str = None
) -> gpd.GeoDataFrame:
    """
    Create a vector shape of the halite data

    Parameters
    ----------
    dat_xr : Xarray dataset of the halite data
    dat_crs : EPSG CRS
    halite : Halite member

    Returns
    -------
    - A (multi)polygon geodataframe of the halite's shape
    """

    if halite:
        shape = (
            dat_xr.sel(halite=halite)["Thickness"]
            .to_dataframe()
            .dropna()
            .reset_index()
        )
    else:
        shape = (
            dat_xr.max(dim="halite")["Thickness"]
            .to_dataframe()
            .dropna()
            .reset_index()
        )
    shape = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries(gpd.points_from_xy(shape.x, shape.y))
        .buffer(100)
        .envelope,
        crs=dat_crs,
    ).dissolve()

    return shape


def zones_of_interest(
    dat_xr: xr.Dataset,
    constraints: dict[str, float],
    dat_crs: int = 23029,
    roof_thickness: float = 80,
    floor_thickness: float = 10,
) -> tuple[gpd.GeoDataFrame, xr.Dataset]:
    """
    Generate a (multi)polygon of the zones of interest by applying thickness
    and depth constraints.

    Parameters
    ----------
    dat_xr : Xarray dataset of the halite data
    constraints : Dictionary containing the following:
        - height: cavern height [m]
        - min_depth: minimum cavern depth [m]
        - max_depth: maximum cavern depth [m]
    dat_crs : EPSG CRS
    roof_thickness : Salt roof thickness [m]
    floor_thickness : Minimum salt floor thickness [m]

    Returns
    -------
    - A (multi)polygon geodataframe of the zones of interest
    - Xarray dataset of the zones of interest
    """

    try:
        zds = dat_xr.where(
            (
                (
                    dat_xr.Thickness
                    >= constraints["height"] + roof_thickness + floor_thickness
                )
                & (
                    dat_xr.TopDepth
                    >= constraints["min_depth"] - roof_thickness
                )
                & (
                    dat_xr.TopDepth
                    <= constraints["max_depth"] - roof_thickness
                )
            ),
            drop=True,
        )
    except KeyError:
        zds = dat_xr.where(
            (
                (
                    dat_xr.Thickness
                    >= constraints["height"] + roof_thickness + floor_thickness
                )
                & (
                    dat_xr.TopDepth
                    >= constraints["min_depth"] - roof_thickness
                )
            ),
            drop=True,
        )

    # zones of interest polygon
    zdf = (
        zds.max(dim="halite")["Thickness"]
        .to_dataframe()
        .dropna()
        .reset_index()
    )
    zdf = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries(gpd.points_from_xy(zdf.x, zdf.y))
        .buffer(100)
        .envelope,
        crs=dat_crs,
    ).dissolve()

    return zdf, zds


def generate_caverns_square_grid(
    dat_extent: gpd.GeoSeries,
    zones_df: gpd.GeoDataFrame,
    diameter: float,
    separation: float = None,
    dat_crs: int = 23029,
) -> gpd.GeoDataFrame:
    """
    Generate salt caverns using a regular square grid within the zones of
    interest.
    Gridding method based on
    https://james-brennan.github.io/posts/fast_gridding_geopandas/.

    Parameters
    ----------
    dat_extent : extent of the data
    zones_df : zones of interest
    diameter : diameter of the cavern [m]
    separation : cavern separation distance [m]
    dat_crs : EPSG CRS

    Returns
    -------
    - A polygon geodataframe of potential caverns
    """

    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds

    if not separation:
        separation = diameter * 4

    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin_, xmax_ + separation, separation):
        for y0 in np.arange(ymin_, ymax_ + separation, separation):
            # bounds
            x1 = x0 - separation
            y1 = y0 + separation
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
    grid_cells = gpd.GeoDataFrame(
        grid_cells, columns=["geometry"], crs=dat_crs
    )

    # generate caverns within the zones of interest
    cavern_df = gpd.sjoin(
        gpd.GeoDataFrame(
            geometry=grid_cells.centroid.buffer(diameter / 2)
        ).drop_duplicates(),
        zones_df,
        predicate="within",
    )

    cavern_df.drop(columns=["index_right"], inplace=True)

    print("Number of potential caverns:", len(cavern_df))

    return cavern_df


def hexgrid_init(
    dat_extent: gpd.GeoSeries, separation: float
) -> tuple[float, list[int], list[int]]:
    """
    Define initial values to be used iteratively to generate the hexagonal
    grid in the `generate_caverns_hexagonal_grid` function

    Parameters
    ----------
    dat_extent : extent of the data
    separation : cavern separation distance [m]

    Returns
    -------
    - a, cols, rows
    """

    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds

    a = np.sin(np.pi / 3)
    cols = np.arange(np.floor(xmin_), np.ceil(xmax_), 3 * separation)
    rows = np.arange(np.floor(ymin_) / a, np.ceil(ymax_) / a, separation)

    return a, cols, rows


def generate_caverns_hexagonal_grid(
    dat_extent: gpd.GeoSeries,
    zones_df: gpd.GeoDataFrame,
    diameter: float,
    separation: float = None,
    dat_crs: int = 23029,
) -> gpd.GeoDataFrame:
    """
    Generate caverns in a regular hexagonal grid as proposed by Williams
    et al. (2022): https://doi.org/10.1016/j.est.2022.105109.
    Hexagonal gridding method based on
    https://sabrinadchan.github.io/data-blog/building-a-hexagonal-cartogram.html.

    Parameters
    ----------
    dat_extent : extent of the data
    zones_df : zones of interest
    diameter : diameter of the cavern [m]
    separation : cavern separation distance [m]
    dat_crs : EPSG CRS

    Returns
    -------
    - A polygon geodataframe of potential caverns
    """

    if not separation:
        separation = diameter * 4

    a, cols, rows = hexgrid_init(dat_extent=dat_extent, separation=separation)

    cavern_df = []
    for x in cols:
        for i, y in enumerate(rows):
            if i % 2 == 0:
                x0 = x
            else:
                x0 = x + 1.5 * separation

            # hexagon vertices
            cavern_df.append(shapely.geometry.Point(x0, y * a))
            cavern_df.append(shapely.geometry.Point(x0 + separation, y * a))
            cavern_df.append(
                shapely.geometry.Point(
                    x0 + (1.5 * separation), (y + separation) * a
                )
            )
            cavern_df.append(
                shapely.geometry.Point(
                    x0 + separation, (y + (2 * separation)) * a
                )
            )
            cavern_df.append(
                shapely.geometry.Point(x0, (y + (2 * separation)) * a)
            )
            cavern_df.append(
                shapely.geometry.Point(
                    x0 - (0.5 * separation), (y + separation) * a
                )
            )
            # hexagon centroid
            cavern_df.append(
                shapely.geometry.Point(
                    x0 + (0.5 * separation), (y + separation) * a
                )
            )

    # generate caverns using hexagon vertices and centroids
    cavern_df = gpd.GeoDataFrame(
        geometry=gpd.GeoDataFrame(geometry=cavern_df, crs=dat_crs)
        .drop_duplicates()
        .buffer(diameter / 2)
    )

    # clip caverns to the zones of interest
    cavern_df = gpd.sjoin(cavern_df, zones_df, predicate="within")

    cavern_df.drop(columns=["index_right"], inplace=True)

    print("Number of potential caverns:", len(cavern_df))

    return cavern_df


def cavern_dataframe(
    dat_zone: xr.Dataset, cavern_df: gpd.GeoDataFrame, dat_crs: int = 23029
) -> gpd.GeoDataFrame:
    """
    Merge halite data for each cavern location

    Parameters
    ----------
    dat_zone : Xarray dataset for the zone of interest
    dat_crs : EPSG CRS
    cavern_df : Geodataframe of caverns within the zone of interest

    Returns
    -------
    - The cavern geodataframe with halite height and depth data for only the
      thickest halite layer at each given point
    """

    # read the halite dataset into a dataframe
    zdf = (
        dat_zone.to_dataframe()[list(dat_zone.data_vars)]
        .dropna()
        .reset_index()
    )

    # generate a square grid for the dataset
    zdf = gpd.GeoDataFrame(
        zdf,
        geometry=gpd.GeoSeries(gpd.points_from_xy(zdf.x, zdf.y))
        .buffer(100)
        .envelope,
        crs=dat_crs,
    )

    # merge with the cavern data
    cavern_df = gpd.sjoin(cavern_df, zdf)

    cavern_df.drop(columns=["index_right"], inplace=True)

    # remove duplicate caverns at each location - keep the thickest layer
    cavern_df = cavern_df.sort_values(
        ["Thickness", "TopDepth"], ascending=False
    ).drop_duplicates(["geometry"])

    return cavern_df


def cavern_volume(height: float, diameter: float, theta: float) -> float:
    """
    Calculate the cavern volume. See Williams et al. (2022), eq. (1) and
    Jannel and Torquet (2022).

    - v_cavern = v_cylinder + 2 * v_cone
    - v_cylinder = pi * r^2 * h_cylinder = pi * r^2 * (h_cavern - 2 * h_cone)
    - v_cone = pi * r^2 * h_cone / 3 = pi * r^2 * (r * tan(theta)) / 3
    - v_cavern = pi * r^2 * (h_cavern - 4 / 3 * r * tan(theta))
    - v_cavern_corrected = v_cavern * scf * (1 - if * insf * bf)

    Parameters
    ----------
    height : Cavern height [m]
    diameter : Cavern diameter [m]
    theta : Cavern roof angle [deg]

    Returns
    -------
    - Corrected cavern volume [m3]
    """

    # calculate ideal cavern volume
    r = diameter / 2
    v_cavern = (
        np.pi * np.square(r) * (height - 4 / 3 * r * np.tan(np.deg2rad(theta)))
    )

    # apply correction factors
    v_cavern = v_cavern * 0.7 * (1 - 0.25 * 0.865 * 1.46)

    return v_cavern


def temperature_cavern_mid_point(
    height: float, depth_top: float, t_0: float = 10, t_delta: float = 27
) -> float:
    """
    Volume of stored hydrogen gas. See Williams et al. (2022), eq. (2).

    Parameters
    ----------
    height : Cavern height [m]
    t_0 : Mean annual surface temperature [deg C]
    t_delta : Geothermal gradient; change in temperature with depth
        [deg C / km]
    depth_top : Cavern top depth [m]

    Returns
    -------
    - Mid point temperature [deg C]
    """

    return t_0 + t_delta / 1000 * (depth_top + height / 2)


def pressure_operating(
    thickness_overburden,
    rho_overburden: float = 2400,
    rho_salt: float = 2200,
    minf: float = 0.3,
    maxf: float = 0.8,
) -> tuple[float, float]:
    """
    Volume of stored hydrogen gas. See Williams et al. (2022), eq. (3) and (4).

    Parameters
    ----------
    thickness_overburden : Overburden thickness / halite top depth [m]
    rho_overburden : Density of the overburden [kg m-3]
    rho_salt : Density of the halite [kg m-3]
    minf : Factor of lithostatic pressure for the minimum operating pressure
    maxf : Factor of lithostatic pressure for the maximum operating pressure

    Returns
    -------
    - Operating pressures [Pa]
    """

    # lithostatic pressure at the casing shoe
    # thickness of the overburden is the same as the depth to top of salt
    # thickness of salt above casing shoe = 50 m
    # acceleration due to gravity = 9.81
    p_casing = (rho_overburden * thickness_overburden + rho_salt * 50) * 9.81

    p_operating_min = minf * p_casing
    p_operating_max = maxf * p_casing

    return p_operating_min, p_operating_max


def density_hydrogen_gas(
    p_operating_min: float,
    p_operating_max: float,
    t_mid_point: float,
    compressibility_factor: float = 1.05,
) -> tuple[float, float]:
    """
    Density of hydrogen at cavern conditions. See Williams et al. (2022),
    s. 3.4.2 and Caglayan et al. (2020), eq. (3).

    http://www.coolprop.org/fluid_properties/fluids/Hydrogen.html

    rho = p * m / (z * r * t)
    where p: pressure [Pa], m: molar mass of hydrogen gas [kg mol-1], z:
    compressibility factor, r: universal gas constant [J K-1 mol-1], t:
    temperature [K]

    Parameters
    ----------
    compressibility_factor : compressibility factor
    p_operating_min : Minimum operating pressure [Pa]
    p_operating_max : Maximum operating pressure [Pa]
    t_mid_point : Mid point temperature [deg C]

    Returns
    -------
    - Hydrogen gas density [kg m-3]
    """

    # rho / p = m / (z * r * t)
    rho_p = (2.01588 / 1000) / (
        compressibility_factor * 8.314 * (t_mid_point + 273.15)
    )

    rho_h2_min = p_operating_min * rho_p
    rho_h2_max = p_operating_max * rho_p

    return rho_h2_min, rho_h2_max


def mass_hydrogen_working(
    rho_h2_min: float, rho_h2_max: float, v_cavern: float
) -> float:
    """
    The working mass of hydrogen in kg. See Williams et al. (2022), eq. (5)
    and (6).

    Parameters
    ----------
    rho_h2_min : Minimum hydrogen density [kg m-3]
    rho_h2_max : Maximum hydrogen density [kg m-3]
    v_cavern : Cavern volume [m3]

    Returns
    -------
    - Working mass of hydrogen [kg]
    """

    # stored mass of hydrogen at min and max operating pressures
    m_min_operating = rho_h2_min * v_cavern
    m_max_operating = rho_h2_max * v_cavern

    return m_max_operating - m_min_operating


def energy_storage_capacity(m_working: float, lhv: float = 119.96) -> float:
    """
    Cavern energy storage capacity in GWh. See Williams et al. (2022), eq. (7).

    Parameters
    ----------
    m_working : Working mass of hydrogen [kg]
    lhv : Lower heating value of hydrogen [MJ kg-1]

    Returns
    -------
    - Energy storage capacity of cavern [GWh]
    """

    return m_working * lhv / 3.6e6
