"""Functions to read and structure data and generate caverns.

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

import glob
import os
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr
from geocube.api.core import make_geocube

# CRS of the Kish Basin dat files
CRS = 23029


def read_dat_file(dat_path):
    """Read XYZ data layers into an Xarray dataset

    Parameters
    ----------
    dat_path : str
        Path to the .dat files

    Returns
    -------
    tuple[xr.Dataset, gpd.GeoSeries]
        Xarray dataset and GeoPandas GeoSeries of the extent
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

    # # find data resolution
    # gdf_xr = gdf[list(gdf.keys())[0]].set_index(["X", "Y"]).to_xarray()
    # resx = gdf_xr["X"][1] - gdf_xr["X"][0]
    # resy = gdf_xr["Y"][1] - gdf_xr["Y"][0]

    # combine dataframes
    gdf = pd.concat(gdf.values())

    # convert dataframe to geodataframe
    gdf["geometry"] = gpd.points_from_xy(gdf.X, gdf.Y, crs=CRS)
    gdf = gpd.GeoDataFrame(gdf)
    gdf.drop(columns=["X", "Y"], inplace=True)

    # convert to Xarray dataset
    xds = make_geocube(
        vector_data=gdf,
        # resolution=(-abs(resy), abs(resx)),
        # align=(abs(resy / 2), abs(resx / 2)),
        resolution=(-200, 200),
        align=(100, 100),
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
        crs=CRS,
    )

    return xds, dat_extent  # zones


def halite_shape(dat_xr, halite=None):
    """Create a vector shape of the halite data

    Parameters
    ----------
    dat_xr : xr.Dataset
        Xarray dataset of the halite data
    halite : str
        Halite member

    Returns
    -------
    gpd.GeoDataFrame
        A (multi)polygon geodataframe of the halite's shape
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
        crs=CRS,
    ).dissolve()

    return shape


def zones_of_interest(
    dat_xr, constraints, roof_thickness=80, floor_thickness=10
):
    """Generate zones of interest by applying thickness and depth constraints.

    Parameters
    ----------
    dat_xr : xr.Dataset
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
    tuple[gpd.GeoDataFrame, xr.Dataset]
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
        crs=CRS,
    ).dissolve()

    return zdf, zds


def generate_caverns_square_grid(
    dat_extent, zones_df, diameter=80, separation=80 * 4
):
    """Generate salt caverns using a regular square grid.

    Parameters
    ----------
    dat_extent : gpd.GeoSeries
        Extent of the data
    zones_df : gpd.GeoDataFrame
        Zones of interest
    diameter : float
        Diameter of the cavern [m]
    separation : float
        Cavern separation distance [m]

    Returns
    -------
    gpd.GeoDataFrame
        A polygon geodataframe of potential caverns in the zone of interest.

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
    grid_cells = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=CRS)

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
    dat_extent : gpd.GeoSeries
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
    dat_extent : gpd.GeoSeries
        Extent of the data
    zones_df : gpd.GeoDataFrame
        Zones of interest
    diameter : float
        Diameter of the cavern [m]
    separation : float
        Cavern separation distance [m]

    Returns
    -------
    gpd.GeoDataFrame
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
        geometry=gpd.GeoDataFrame(geometry=cavern_df, crs=CRS)
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
    dat_zone : xr.Dataset
        Xarray dataset for the zone of interest
    cavern_df : gpd.GeoDataFrame
        Geodataframe of caverns within the zone of interest

    Returns
    -------
    gpd.GeoDataFrame
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
        crs=CRS,
    )

    # merge with the cavern data
    cavern_df = gpd.sjoin(cavern_df, zdf)

    cavern_df.drop(columns=["index_right"], inplace=True)

    # remove duplicate caverns at each location - keep the thickest layer
    cavern_df = cavern_df.sort_values(
        ["Thickness", "TopDepth"], ascending=False
    ).drop_duplicates(["geometry"])

    return cavern_df


def read_shapefile_from_zip(data_path, endswith=".shp"):
    """Read the Shapefile layer from a Zip file.

    Parameters
    ----------
    data_path : str
        Path to the exploration well Zip file
    endswith : str
        What the Shapefile's filename ends with

    Returns
    -------
    gpd.GeoDataFrame
        Geodataframe of the Shapefile's data
    """
    data_shp = gpd.read_file(
        os.path.join(
            f"zip://{data_path}!"
            + [
                x
                for x in ZipFile(data_path).namelist()
                if x.endswith(endswith)
            ][0]
        )
    )

    return data_shp


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
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    500 m buffer - suggested in the draft OREDP II p. 108 [#DECC23]_.
    """
    wells = read_shapefile_from_zip(data_path=os.path.join(data_path))

    wells = wells[wells["AREA"].str.contains("Kish")].to_crs(CRS)

    wells_b = gpd.GeoDataFrame(geometry=wells.buffer(buffer))

    return wells, wells_b


def constraint_wind_farm(data_path, dat_extent):
    """Read data for wind farms.

    Parameters
    ----------
    data_path : str
        Path to the Zip file
    dat_extent : gpd.GeoSeries
        Extent of the data

    Returns
    -------
    gpd.GeoDataFrame
        Geodataframes of the dataset

    Notes
    -----
    The shapes are used as is without a buffer - suggested for renewable
    energy test site areas in the draft OREDP II p. 109 [#DECC23]_.
    """
    wind_farms = read_shapefile_from_zip(data_path=os.path.join(data_path))

    # keep only features near Kish Basin
    wind_farms = (
        wind_farms.to_crs(CRS)
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
    dat_extent : gpd.GeoSeries
        Extent of the data
    buffer : float
        Buffer [m]

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    1 NM (nautical mile) buffer - suggested in the draft OREDP II p. 108
    [#DECC23]_. 1 NM is equivalent to 1,852 m [#NIST16]_.
    """
    shipping = read_shapefile_from_zip(data_path=os.path.join(data_path))

    # keep only features near Kish Basin
    shipping = (
        shipping.to_crs(CRS)
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
    dat_extent : gpd.GeoSeries
        Extent of the data
    buffer : float
        Buffer [m]

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    Archaeological Exclusion Zones recommendation - 100 m buffer [#RE21]_.
    """
    shipwrecks = read_shapefile_from_zip(data_path=os.path.join(data_path))

    # keep only features near Kish Basin
    shipwrecks = (
        shipwrecks.to_crs(CRS)
        .sjoin(gpd.GeoDataFrame(geometry=dat_extent.buffer(3000)))
        .reset_index()
    )

    shipwrecks.drop(columns=["index_right"], inplace=True)

    # Archaeological Exclusion Zones recommendation
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
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        Geodataframes of the dataset and buffer

    Notes
    -----
    750 m buffer - suggested in the draft OREDP II p. 109-111 [#DECC23]_.
    """
    cables = gpd.read_file(os.path.join(data_path))

    cables = cables.to_crs(CRS)

    # remove point features
    cables = cables.drop(range(3)).reset_index(drop=True)

    cables_b = gpd.GeoDataFrame(geometry=cables.buffer(buffer)).dissolve()

    return cables, cables_b


def constraint_halite_edge(dat_xr, buffer=80 * 3):
    """The edge of each halite member as a constraint.

    Parameters
    ----------
    dat_xr : xr.Dataset
        Xarray dataset of the halite data
    buffer : float
        Buffer [m]

    Returns
    -------
    dict[str, gpd.GeoDataFrame]
        Dictionary of GeoPandas geodataframes of the halite edge constraint
        for each halite member

    Notes
    -----
    Set the buffer to 3 times the cavern diameter, i.e. the pillar width.
    """

    buffer_edge = {}
    for halite in dat_xr.halite.values:
        buffer_edge[halite] = halite_shape(dat_xr=dat_xr, halite=halite)
        buffer_edge[halite] = gpd.GeoDataFrame(
            geometry=buffer_edge[halite].boundary.buffer(buffer)
        ).overlay(buffer_edge[halite], how="intersection")

    return buffer_edge


def generate_caverns_with_constraints(
    zones_gdf, zones_ds, dat_extent, exclusions, diameter=80, separation=80 * 4
):
    """Add constraints to cavern configuration.

    Parameters
    ----------
    zones_gdf : gpd.GeoDataFrame
        GeoPandas dataframe of zones of interest
    zones_ds : xr.Dataset
        Xarray dataset of zones of interest
    dat_extent : gpd.GeoSeries
        Extent of the data
    exclusions : dict[str, gpd.GeoDataFrame]
        A dictionary of exclusions data. If any of the following keys do not
        exist in the dictionary, the exclusion will be skipped:
        `"edge"`: halite edge, dict[str, gpd.GeoDataFrame];
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
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
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
    try:
        cavern_df = cavern_df.overlay(
            gpd.sjoin(
                cavern_df, exclusions["shipping"], predicate="intersects"
            ),
            how="difference",
        )
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")

    print("-" * 60)
    print("Exclude subsea cables...")
    try:
        cavern_df = cavern_df.overlay(
            gpd.sjoin(cavern_df, exclusions["cables"], predicate="intersects"),
            how="difference",
        )
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")

    print("-" * 60)
    print("Exclude wind farms...")
    try:
        cavern_df = cavern_df.overlay(
            gpd.sjoin(
                cavern_df, exclusions["wind_farms"], predicate="intersects"
            ),
            how="difference",
        )
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")

    print("-" * 60)
    print("Exclude exploration wells...")
    try:
        cavern_df = cavern_df.overlay(
            gpd.sjoin(cavern_df, exclusions["wells"], predicate="intersects"),
            how="difference",
        )
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")

    print("-" * 60)
    print("Exclude shipwrecks...")
    try:
        cavern_df = cavern_df.overlay(
            gpd.sjoin(
                cavern_df, exclusions["shipwrecks"], predicate="intersects"
            ),
            how="difference",
        )
        print(f"Number of potential caverns: {len(cavern_df):,}")
        pct = (len(cavern_all) - len(cavern_df)) / len(cavern_all) * 100
        print(f"Caverns excluded: {pct:.2f}%")
    except KeyError:
        print("No data specified!")

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
    cavern_df : gpd.GeoDataFrame
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
    gpd.GeoDataFrame
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
