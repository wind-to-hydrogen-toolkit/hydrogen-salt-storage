"""Functions to structure data and generate caverns with constraints.

.. rubric:: References
.. [#Brennan20] Brennan, J. (2020). ‘Fast and easy gridding of point data with
    geopandas’, 16 March. Available at:
    https://james-brennan.github.io/posts/fast_gridding_geopandas/
    (Accessed: 1 January 2024).
.. [#Williams22] Williams, J. D. O., Williamson, J. P., Parkes, D., Evans, D.
    J., Kirk, K. L., Sunny, N., Hough, E., Vosper, H., and Akhurst, M. C.
    (2022). ‘Does the United Kingdom have sufficient geological storage
    capacity to support a hydrogen economy? Estimating the salt cavern storage
    potential of bedded halite formations’, Journal of Energy Storage, 53, p.
    105109. https://doi.org/10.1016/j.est.2022.105109.
.. [#Chan19] Chan, S. (2019). ‘Left Join – Building a hexagonal cartogram
    with Python’, 7 September. Available at:
    https://sabrinadchan.github.io/data-blog/building-a-hexagonal-cartogram.html
    (Accessed: 1 January 2024).
.. [#DECC23] Department of the Environment, Climate and Communications (2023).
    Draft Offshore Renewable Energy Development Plan II (OREDP II). Government
    of Ireland. Available at:
    https://www.gov.ie/en/publication/71e36-offshore-renewable-energy-development-plan-ii-oredp-ii/
    (Accessed: 1 October 2023).
.. [#NIST16] National Institute of Standards and Technology (2016). ‘Units
    Outside the SI’, in NIST Guide to the SI. (Special Publication, 811).
    Available at:
    https://www.nist.gov/pml/special-publication-811/nist-guide-si-chapter-5-units-outside-si
    (Accessed: 30 November 2023).
.. [#RE21] Renewables Ireland (2021). Dublin Array Offshore Wind Farm:
    Supporting Information Report - Annex D: Marine Archaeology Assessment.
    Available at:
    https://assets.gov.ie/202763/9160dd51-b119-44cf-a224-8b5302347e7d.pdf
    (Accessed: 10 November 2023).
"""

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from src import read_data as rd


def zones_of_interest(
    dat_xr, constraints, roof_thickness=80, floor_thickness=10
):
    """Generate zones of interest by applying thickness and depth constraints.

    Parameters
    ----------
    dat_xr : xarray.Dataset
        Xarray dataset of the halite data
    constraints : dict[str, float]
        Dictionary containing the following:
        `"height"`: cavern height [m];
        `"min_depth"`: minimum cavern depth [m];
        `"max_depth"`: Maximum cavern depth [m]
    roof_thickness : float
        Salt roof thickness [m]
    floor_thickness : float
        Minimum salt floor thickness [m]

    Returns
    -------
    tuple[geopandas.GeoDataFrame, xarray.Dataset]
        A (multi)polygon geodataframe of the zones of interest and an Xarray
        dataset of the zones of interest
    """
    zds = dat_xr.where(
        (
            (
                dat_xr.Thickness
                >= constraints["height"] + roof_thickness + floor_thickness
            )
            & (dat_xr.TopDepth >= constraints["min_depth"] - roof_thickness)
            & (dat_xr.TopDepth <= constraints["max_depth"] - roof_thickness)
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
        crs=rd.CRS,
    ).dissolve()

    return zdf, zds


def generate_caverns_square_grid(
    dat_extent, zones_df, diameter=80, separation=80 * 4
):
    """Generate salt caverns using a regular square grid.

    Parameters
    ----------
    dat_extent : geopandas.GeoSeries
        Extent of the data
    zones_df : geopandas.GeoDataFrame
        Zones of interest
    diameter : float
        Diameter of the cavern [m]
    separation : float
        Cavern separation distance [m]

    Returns
    -------
    geopandas.GeoDataFrame
        A polygon geodataframe of potential caverns in the zone of interest

    Notes
    -----
    Gridding method based on [#Brennan20]_.
    """
    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds

    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin_, xmax_ + separation, separation):
        for y0 in np.arange(ymin_, ymax_ + separation, separation):
            # bounds
            x1 = x0 - separation
            y1 = y0 + separation
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
    grid_cells = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=rd.CRS)

    # generate caverns within the zones of interest
    cavern_df = gpd.sjoin(
        gpd.GeoDataFrame(
            geometry=grid_cells.centroid.buffer(diameter / 2)
        ).drop_duplicates(),
        zones_df,
        predicate="within",
    )

    cavern_df.drop(columns=["index_right"], inplace=True)

    print(f"Number of potential caverns: {len(cavern_df):,}")

    return cavern_df


def hexgrid_init(dat_extent, separation):
    """Initialise variables for `generate_caverns_hexagonal_grid`.

    Parameters
    ----------
    dat_extent : geopandas.GeoSeries
        Extent of the data
    separation : float
        Cavern separation distance [m]

    Returns
    -------
    tuple[float, list[int], list[int]]
        a, cols, and rows
    """
    xmin_, ymin_, xmax_, ymax_ = dat_extent.total_bounds
    a = np.sin(np.pi / 3)
    cols = np.arange(np.floor(xmin_), np.ceil(xmax_), 3 * separation)
    rows = np.arange(np.floor(ymin_) / a, np.ceil(ymax_) / a, separation)
    return a, cols, rows


def generate_caverns_hexagonal_grid(
    dat_extent, zones_df, diameter=80, separation=80 * 4
):
    """Generate caverns in a regular hexagonal grid.

    Parameters
    ----------
    dat_extent : geopandas.GeoSeries
        Extent of the data
    zones_df : geopandas.GeoDataFrame
        Zones of interest
    diameter : float
        Diameter of the cavern [m]
    separation : float
        Cavern separation distance [m]

    Returns
    -------
    geopandas.GeoDataFrame
        A polygon geodataframe of potential caverns

    Notes
    -----
    The close-packed hexagonal grid configuration was proposed by
    [#Williams22]_; this configuration provides around 15% more caverns
    compared to a square grid configuration. Hexagonal gridding method based
    on [#Chan19]_.
    """
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
        geometry=gpd.GeoDataFrame(geometry=cavern_df, crs=rd.CRS)
        .drop_duplicates()
        .buffer(diameter / 2)
    )

    # clip caverns to the zones of interest
    cavern_df = gpd.sjoin(cavern_df, zones_df, predicate="within")

    cavern_df.drop(columns=["index_right"], inplace=True)

    print(f"Number of potential caverns: {len(cavern_df):,}")

    return cavern_df


def cavern_dataframe(dat_zone, cavern_df):
    """Merge halite data for each cavern location and create a dataframe.

    Parameters
    ----------
    dat_zone : xarray.Dataset
        Xarray dataset for the zone of interest
    cavern_df : geopandas.GeoDataFrame
        Geodataframe of caverns within the zone of interest

    Returns
    -------
    geopandas.GeoDataFrame
        The cavern geodataframe with halite height and depth data for only the
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
        crs=rd.CRS,
    )

    # merge with the cavern data
    cavern_df = gpd.sjoin(cavern_df, zdf)

    cavern_df.drop(columns=["index_right"], inplace=True)

    # remove duplicate caverns at each location - keep the thickest layer
    cavern_df = cavern_df.sort_values(
        ["Thickness", "TopDepth"], ascending=False
    ).drop_duplicates(["geometry"])

    return cavern_df


def constraint_exploration_well(data_path, buffer=500):
    """Read exploration well data and generate constraint.

    Parameters
    ----------
    data_path : str
        Path to the Zip file
    buffer : float
        Buffer [m]

    Returns
    -------
    tuple[geopandas.GeoDataFrame, geopandas.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    500 m buffer - suggested in the draft OREDP II p. 108 [#DECC23]_.
    """
    wells = rd.read_shapefile_from_zip(data_path=os.path.join(data_path))
    wells = wells[wells["AREA"].str.contains("Kish")].to_crs(rd.CRS)
    wells_b = gpd.GeoDataFrame(geometry=wells.buffer(buffer))
    return wells, wells_b


def constraint_wind_farm(data_path, dat_extent):
    """Read data for wind farms.

    Parameters
    ----------
    data_path : str
        Path to the Zip file
    dat_extent : geopandas.GeoSeries
        Extent of the data

    Returns
    -------
    geopandas.GeoDataFrame
        Geodataframe of the dataset

    Notes
    -----
    The shapes are used as is without a buffer - suggested for renewable
    energy test site areas in the draft OREDP II p. 109 [#DECC23]_.
    """
    wind_farms = rd.read_shapefile_from_zip(data_path=os.path.join(data_path))
    # keep only features near Kish Basin
    wind_farms = (
        wind_farms.to_crs(rd.CRS)
        .sjoin(gpd.GeoDataFrame(geometry=dat_extent.buffer(3000)))
        .reset_index()
        .sort_values("Name")
    )
    wind_farms.drop(columns=["index_right"], inplace=True)
    return wind_farms


def constraint_shipping_routes(data_path, dat_extent, buffer=1852):
    """Read frequent shipping route data and generate constraint.

    Parameters
    ----------
    data_path : str
        Path to the Zip file
    dat_extent : geopandas.GeoSeries
        Extent of the data
    buffer : float
        Buffer [m]

    Returns
    -------
    tuple[geopandas.GeoDataFrame, geopandas.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    1 NM (nautical mile) buffer - suggested in the draft OREDP II p. 108
    [#DECC23]_. 1 NM is equivalent to 1,852 m [#NIST16]_.
    """
    shipping = rd.read_shapefile_from_zip(data_path=os.path.join(data_path))
    # keep only features near Kish Basin
    shipping = (
        shipping.to_crs(rd.CRS)
        .sjoin(gpd.GeoDataFrame(geometry=dat_extent.buffer(3000)))
        .reset_index()
    )
    shipping.drop(columns=["index_right"], inplace=True)
    shipping_b = gpd.GeoDataFrame(geometry=shipping.buffer(buffer)).dissolve()
    return shipping, shipping_b


def constraint_shipwrecks(data_path, dat_extent, buffer=100):
    """Read shipwreck data and generate constraint.

    Parameters
    ----------
    data_path : str
        Path to the Zip file
    dat_extent : geopandas.GeoSeries
        Extent of the data
    buffer : float
        Buffer [m]

    Returns
    -------
    tuple[geopandas.GeoDataFrame, geopandas.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    Archaeological Exclusion Zones recommendation - 100 m buffer [#RE21]_.
    """
    shipwrecks = rd.read_shapefile_from_zip(data_path=os.path.join(data_path))
    # keep only features near Kish Basin
    shipwrecks = (
        shipwrecks.to_crs(rd.CRS)
        .sjoin(gpd.GeoDataFrame(geometry=dat_extent.buffer(3000)))
        .reset_index()
    )
    shipwrecks.drop(columns=["index_right"], inplace=True)
    shipwrecks_b = gpd.GeoDataFrame(geometry=shipwrecks.buffer(buffer))
    return shipwrecks, shipwrecks_b


def constraint_subsea_cables(data_path, buffer=750):
    """Read subsea cable data and generate constraint.

    Parameters
    ----------
    data_path : str
        Path to the GPKG file
    buffer : float
        Buffer [m]

    Returns
    -------
    tuple[geopandas.GeoDataFrame, geopandas.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    750 m buffer - suggested in the draft OREDP II p. 109-111 [#DECC23]_.
    """
    cables = gpd.read_file(os.path.join(data_path))
    cables = cables.to_crs(rd.CRS)
    # remove point features
    cables = cables.drop(range(3)).reset_index(drop=True)
    cables_b = gpd.GeoDataFrame(geometry=cables.buffer(buffer)).dissolve()
    return cables, cables_b


def constraint_halite_edge(dat_xr, buffer=80 * 3):
    """The edge of each halite member as a constraint.

    Parameters
    ----------
    dat_xr : xarray.Dataset
        Xarray dataset of the halite data
    buffer : float
        Buffer [m]

    Returns
    -------
    dict[str, geopandas.GeoDataFrame]
        Dictionary of GeoPandas geodataframes of the halite edge constraint
        for each halite member

    Notes
    -----
    Set the buffer to 3 times the cavern diameter, i.e. the pillar width.
    """
    buffer_edge = {}
    for halite in dat_xr.halite.values:
        buffer_edge[halite] = rd.halite_shape(dat_xr=dat_xr, halite=halite)
        buffer_edge[halite] = gpd.GeoDataFrame(
            geometry=buffer_edge[halite].boundary.buffer(buffer)
        ).overlay(buffer_edge[halite], how="intersection")
    return buffer_edge


def exclude_constraint(cavern_df, cavern_all, exclusions, key):
    """Exclude constraint by their dictionary key.

    Parameters
    ----------
    cavern_df : geopandas.GeoDataFrame
        Dataframe of available caverns
    cavern_all : geopandas.GeoDataFrame
        Dataframe of all caverns, i.e. available and excluded
    exclusions : dict[str, geopandas.GeoDataFrame]
        Dictionary of exclusions data
    key : str
        Key for the constraint in the dictionary; one of the following:
        `"shipping"`: frequent shipping routes;
        `"cables"`: subsea cables;
        `"wind_farms"`: offshore wind farms;
        `"wells"`: exporation wells;
        `"shipwrecks"`: shipwrecks

    Returns
    -------
    geopandas.GeoDataFrame
        Dataframe of available caverns
    """
    try:
        cavern_df = cavern_df.overlay(
            gpd.sjoin(cavern_df, exclusions[key], predicate="intersects"),
            how="difference",
        )
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")
    print("-" * 60)
    return cavern_df


def generate_caverns_with_constraints(
    zones_gdf, zones_ds, dat_extent, exclusions, diameter=80, separation=80 * 4
):
    """Add constraints to cavern configuration.

    Parameters
    ----------
    zones_gdf : geopandas.GeoDataFrame
        GeoPandas dataframe of zones of interest
    zones_ds : xarray.Dataset
        Xarray dataset of zones of interest
    dat_extent : geopandas.GeoSeries
        Extent of the data
    exclusions : dict[str, geopandas.GeoDataFrame]
        Dictionary of exclusions data; if any of the following keys do not
        exist in the dictionary, the exclusion will be skipped:
        `"edge"`: halite edge, dict[str, geopandas.GeoDataFrame];
        `"shipping"`: frequent shipping routes;
        `"cables"`: subsea cables;
        `"wind_farms"`: offshore wind farms;
        `"wells"`: exporation wells;
        `"shipwrecks"`: shipwrecks
    diameter : float
        Diameter of the cavern [m]
    separation : float
        Cavern separation distance [m]

    Returns
    -------
    tuple[geopandas.GeoDataFrame, geopandas.GeoDataFrame]
        Dataframe of available and excluded caverns
    """
    print("Without constraints...")
    cavern_df = generate_caverns_hexagonal_grid(
        dat_extent=dat_extent,
        zones_df=zones_gdf,
        diameter=diameter,
        separation=separation,
    )
    cavern_df = cavern_dataframe(dat_zone=zones_ds, cavern_df=cavern_df)
    # keep original
    cavern_all = cavern_df.copy()
    print("-" * 60)

    print("Exclude salt formation edges...")
    try:
        cavern_dict = {}
        for h in cavern_df["halite"].unique():
            cavern_dict[h] = cavern_df[cavern_df["halite"] == h]
            cavern_dict[h] = cavern_dict[h].overlay(
                gpd.sjoin(
                    cavern_dict[h],
                    exclusions["edge"][h],
                    predicate="intersects",
                ),
                how="difference",
            )
        cavern_df = pd.concat(cavern_dict.values())
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")
    print("-" * 60)

    print("Exclude frequent shipping routes...")
    cavern_df = exclude_constraint(
        cavern_df=cavern_df,
        cavern_all=cavern_all,
        exclusions=exclusions,
        key="shipping",
    )

    print("Exclude subsea cables...")
    cavern_df = exclude_constraint(
        cavern_df=cavern_df,
        cavern_all=cavern_all,
        exclusions=exclusions,
        key="cables",
    )

    print("Exclude wind farms...")
    cavern_df = exclude_constraint(
        cavern_df=cavern_df,
        cavern_all=cavern_all,
        exclusions=exclusions,
        key="wind_farms",
    )

    print("Exclude exploration wells...")
    cavern_df = exclude_constraint(
        cavern_df=cavern_df,
        cavern_all=cavern_all,
        exclusions=exclusions,
        key="wells",
    )

    print("Exclude shipwrecks...")
    cavern_df = exclude_constraint(
        cavern_df=cavern_df,
        cavern_all=cavern_all,
        exclusions=exclusions,
        key="shipwrecks",
    )

    # get excluded caverns
    caverns_excl = cavern_all.overlay(
        gpd.sjoin(cavern_all, cavern_df, predicate="intersects"),
        how="difference",
    )

    return cavern_df, caverns_excl


def label_caverns(
    cavern_df, heights, depths, roof_thickness=80, floor_thickness=10
):
    """Label cavern dataframe by height and depth.

    Parameters
    ----------
    cavern_df : geopandas.GeoDataFrame
        Dataframe of potential caverns
    heights : list[float]
        List of fixed caverns heights [m] for labelling
    depths : dict[str, float]
        Dictionary of cavern top depth ranges [m] for labelling
    roof_thickness : float
        Salt roof thickness [m]
    floor_thickness : float
        Minimum salt floor thickness [m]

    Returns
    -------
    geopandas.GeoDataFrame
        A dataframe of potential caverns labelled by cavern height and top
        depth ranges
    """
    # label caverns by height
    if len(heights) == 1:
        cavern_df["height"] = heights[0]
    else:
        conditions = [
            (
                cavern_df["Thickness"]
                < (heights[1] + roof_thickness + floor_thickness)
            ),
            (
                cavern_df["Thickness"]
                >= (heights[1] + roof_thickness + floor_thickness)
            )
            & (
                cavern_df["Thickness"]
                < (heights[2] + roof_thickness + floor_thickness)
            ),
            (
                cavern_df["Thickness"]
                >= (heights[2] + roof_thickness + floor_thickness)
            ),
        ]
        choices = [str(x) for x in heights]
        cavern_df["height"] = np.select(conditions, choices)

    # label caverns by depth
    conditions = [
        (cavern_df["TopDepth"] < (depths["min_opt"] - roof_thickness)),
        (cavern_df["TopDepth"] >= (depths["min_opt"] - roof_thickness))
        & (cavern_df["TopDepth"] <= (depths["max_opt"] - roof_thickness)),
        (cavern_df["TopDepth"] > (depths["max_opt"] - roof_thickness)),
    ]
    choices = [
        f"{depths['min']:,} - {depths['min_opt']:,}",
        f"{depths['min_opt']:,} - {depths['max_opt']:,}",
        f"{depths['max_opt']:,} - {depths['max']:,}",
    ]
    cavern_df["depth"] = np.select(conditions, choices)

    # create columns for the cavern heights and top depths
    cavern_df["cavern_height"] = cavern_df["height"].astype(float)
    cavern_df["cavern_depth"] = cavern_df["TopDepth"] + roof_thickness

    return cavern_df


def read_weibull_data(data_path_weibull, data_path_wind_farms, dat_extent):
    """Extract mean, max, and min Weibull parameters of wind speeds.

    Parameters
    ----------
    data_path_weibull : str
        Path to the Weibull parameter data Zip file
    data_path_wind_farms : str
        Path to the wind farm data Zip file
    dat_extent : geopandas.GeoSeries
        Extent of the Kish Basin data

    Returns
    -------
    pandas.DataFrame
        Dataframe of k and C values for each wind farm

    Notes
    -----
    Data extracted for each wind farm in the area of interest, i.e. Kish
    Basin: Codling, Dublin Array, and NISA.
    """
    weibull_df = {}

    for w in ["c", "k"]:
        # read Weibull parameter data
        weibull_df[w] = rd.read_shapefile_from_zip(
            data_path=data_path_weibull, endswith=f"{w}_ITM.shp"
        )

        # read wind farm data
        wind_farms = constraint_wind_farm(
            data_path=data_path_wind_farms,
            dat_extent=dat_extent,
        )

        # for combining Codling wind farm polygons
        wind_farms["Name_"] = wind_farms["Name"].str.split(expand=True)[0]

        # convert CRS and keep areas intersecting with wind farms
        weibull_df[w] = (
            weibull_df[w]
            .to_crs(rd.CRS)
            .overlay(wind_farms, how="intersection")
        )

        # rename column
        weibull_df[w].rename(columns={"Value": w}, inplace=True)

        # average c and k over wind farms
        weibull_df[w] = wind_farms.dissolve(by="Name_").merge(
            weibull_df[w].dissolve(
                by="Name_", aggfunc={w: ["min", "max", "mean"]}
            ),
            on="Name_",
        )

        # keep only relevant columns
        weibull_df[w] = weibull_df[w][
            ["Name", (w, "min"), (w, "max"), (w, "mean")]
        ]

        # reset index
        weibull_df[w] = weibull_df[w].reset_index(drop=True)

    # merge
    weibull_df = pd.merge(weibull_df["c"], weibull_df["k"], on="Name")

    return weibull_df
