"""Functions to download and read data.

"""

import glob
import os
from datetime import datetime, timezone
from zipfile import ZipFile

import geopandas as gpd
import pandas as pd
import pooch
import rasterio as rio
import shapely
import xarray as xr
from geocube.api.core import make_geocube

# CRS of the Kish Basin dat files
CRS = 23029


def download_data(url, data_dir, file_name, known_hash=None):
    """Download data and store it in the specified directory using Pooch.

    Parameters
    ----------
    url : str
        URL from which the data will be downloaded
    data_dir : str
        Directory to store the downloaded data
    file_name : str
        Name of the downloaded data file with its extension (not full path)
    known_hash : str
        SHA256 hash of downloaded file

    Notes
    -----
    This only downloads data if necessary, i.e. if the data file does not
    already exist in the directory.
    """
    data_file = os.path.join(data_dir, file_name)
    if not os.path.isfile(data_file):
        os.makedirs(data_dir, exist_ok=True)
        pooch.retrieve(
            url=url, known_hash=known_hash, fname=file_name, path=data_dir
        )
        print(f"Data downloaded on: {datetime.now(tz=timezone.utc)}")
        with open(f"{data_file}.txt", "w", encoding="utf-8") as outfile:
            outfile.write(
                f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
                f"Download URL: {url}\n"
                f"SHA256 hash: {pooch.file_hash(data_file)}\n"
            )
    else:
        print(f"Data '{file_name}' already exists in '{data_dir}'.")
        with open(f"{data_file}.txt", encoding="utf-8") as f:
            print(f.read())


def kish_basin_extent(dat_path):
    """Read the extent of the Kish Basin data.

    Parameters
    ----------
    dat_path : str
        Path to the Kish Basin data files

    Returns
    -------
    geopandas.GeoSeries
        GeoPandas geoseries of the extent polygon
    """
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
    return dat_extent


def read_dat_file(dat_path):
    """Read XYZ halite data layers into an Xarray dataset.

    Parameters
    ----------
    dat_path : str
        Path to the DAT (ASCII XYZ) files

    Returns
    -------
    tuple[xarray.Dataset, geopandas.GeoSeries]
        Xarray dataset of the XYZ data and GeoPandas geoseries of the extent
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
        # unit = d.split(" ")[-1]
        zvar = d.split("Halite ")[-1].split(" XYZ")[0]
        xds_[d] = (
            xds.sel(data=d)
            .assign_coords(halite=halite_member)
            .expand_dims(dim="halite")
            .drop_vars("data")
        )
        xds_[d] = xds_[d].rename({"Z": zvar.replace(" ", "")})
        # xds_[d][zvar.replace(" ", "")] = xds_[d][
        #     zvar.replace(" ", "")
        # ].assign_attrs(units=unit, long_name=zvar)

    xds = xr.combine_by_coords(xds_.values(), combine_attrs="override")

    # # keep only points corresponding to zones of interest in the dataframe
    # zones = gdf.loc[gdf["data"].str.contains("Zone")]

    # # create zones of interest polygon
    # zones = gpd.GeoDataFrame(geometry=zones.buffer(100).envelope).dissolve()

    # create extent polygon
    dat_extent = kish_basin_extent(dat_path=dat_path)

    return xds, dat_extent  # zones


def kish_basin_data_depth_adjusted(dat_path, bathymetry_path):
    """Read halite data with adjusted depths and its extent.

    Parameters
    ----------
    dat_path : str
        Path to the DAT (ASCII XYZ) files
    bathymetry_path : str
        Path to the bathymetry netCDF file

    Returns
    -------
    tuple[xarray.Dataset, geopandas.GeoSeries]
        Xarray dataset of the halite data and GeoPandas geoseries of the extent

    Notes
    -----
    A hacky workaround was used to prevent multiple grid mappings in the
    resulting Xarray dataset by first assigning a variable of the desired
    mapping and then multiplying it by zero prior to involving variables which
    originally had a different mapping.
    """
    bath = xr.open_dataset(
        os.path.join(bathymetry_path, "D4_2022.nc"), decode_coords="all"
    )
    bath = bath.chunk({"lat": 1000, "lon": 1000, "cdi_index_count": 1000})
    xds, dat_extent = read_dat_file(dat_path=dat_path)
    bath = bath.rename({"lon": "x", "lat": "y"}).rio.reproject_match(
        xds, resampling=rio.enums.Resampling.bilinear
    )
    xds = xds.assign(TopDepthSeabed=xds["TopDepth"] + bath["elevation"])
    xds = xds.assign(BaseDepthSeabed=xds["BaseDepth"] + bath["elevation"])
    # hacky way to prevent multiple grid mappings
    xds = xds.assign(Bathymetry=xds["TopDepth"] * 0 + bath["elevation"])
    # for v in ["TopDepth", "BaseDepth"]:
    #     xds[f"{v}Seabed"].attrs = xds[v].attrs
    #     xds[f"{v}Seabed"].attrs["long_name"] = (
    #         xds[f"{v}Seabed"].attrs["long_name"] + " relative to the "
    #         "sea-floor"
    #     )
    # xds["Bathymetry"].attrs = bath["elevation"].attrs
    return xds, dat_extent


def bathymetry_layer(dat_extent, bathymetry_path):
    """Bathymetry layer for plotting.

    Parameters
    ----------
    dat_extent : geopandas.GeoSeries
        Extent of the Kish Basin data
    bathymetry_path : str
        Path to the bathymetry netCDF file

    Returns
    -------
    xarray.Dataset
        Xarray dataset of the halite data
    """
    bath = xr.open_dataset(
        os.path.join(bathymetry_path, "D4_2022.nc"), decode_coords="all"
    )
    bath = bath.chunk({"lat": 1000, "lon": 1000, "cdi_index_count": 1000})
    bath = bath.rio.reproject(CRS).rio.clip(dat_extent.buffer(3000))
    return bath


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
    geopandas.GeoDataFrame
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


def halite_shape(dat_xr, halite=None):
    """Create a vector shape of the halite data.

    Parameters
    ----------
    dat_xr : xarray.Dataset
        Xarray dataset of the halite data
    halite : str
        Halite member

    Returns
    -------
    geopandas.GeoDataFrame
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
