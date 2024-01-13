#!/usr/bin/env python
# coding: utf-8

# # Wind farm optimisation

import os

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import Point

from h2ss import capacity as cap
from h2ss import data as rd
from h2ss import functions as fns
from h2ss import optimisation as opt

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# ## Halite data

ds, extent = rd.kish_basin_data_depth_adjusted(
    dat_path=os.path.join("data", "kish-basin"),
    bathymetry_path=os.path.join("data", "bathymetry"),
)

# ## Constraints

# exploration wells
_, wells_b = fns.constraint_exploration_well(
    data_path=os.path.join(
        "data",
        "exploration-wells",
        "Exploration_Wells_Irish_Offshore.shapezip.zip",
    )
)

# wind farms
wind_farms = fns.constraint_wind_farm(
    data_path=os.path.join(
        "data", "wind-farms", "wind-farms-foreshore-process.zip"
    ),
    dat_extent=extent,
)

# frequent shipping routes
_, shipping_b = fns.constraint_shipping_routes(
    data_path=os.path.join(
        "data", "shipping", "shipping_frequently_used_routes.zip"
    ),
    dat_extent=extent,
)

# shipwrecks
_, shipwrecks_b = fns.constraint_shipwrecks(
    data_path=os.path.join(
        "data", "shipwrecks", "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"
    ),
    dat_extent=extent,
)

# subsea cables
_, cables_b = fns.constraint_subsea_cables(
    data_path=os.path.join("data", "subsea-cables", "KIS-ORCA.gpkg")
)

# distance from salt formation edge
edge_buffer = fns.constraint_halite_edge(dat_xr=ds)

# ## Zones of interest

# height = 155 m, 1,000 m <= depth <= 1,500 m, diameter = 80 m,
# separation = 320 m
zones, zds = fns.zones_of_interest(
    dat_xr=ds,
    constraints={"height": 155, "min_depth": 1000, "max_depth": 1500},
)

# ## Exclusions

caverns, _ = fns.generate_caverns_with_constraints(
    zones_gdf=zones,
    zones_ds=zds,
    dat_extent=extent,
    exclusions={
        "wells": wells_b,
        "wind_farms": wind_farms,
        "shipwrecks": shipwrecks_b,
        "shipping": shipping_b,
        "cables": cables_b,
        "edge": edge_buffer,
    },
)

# label caverns by height and depth
caverns = fns.label_caverns(
    cavern_df=caverns,
    heights=[155],
    depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
)

# ## Capacity [GWh]

# calculate volumes and capacities
caverns = cap.calculate_capacity_dataframe(cavern_df=caverns)

# total capacity
caverns[["capacity"]].sum().iloc[0]

# total working mass
caverns[["working_mass"]].sum().iloc[0]

# ## Power curve [MW] and Weibull wind speed distribution

# extract data for wind farms at 150 m
data = fns.read_weibull_data(
    data_path_weibull=os.path.join(
        "data", "weibull-parameters-wind-speeds", "Weibull_150m_params_ITM.zip"
    ),
    data_path_wind_farms=os.path.join(
        "data", "wind-farms", "wind-farms-foreshore-process.zip"
    ),
    dat_extent=extent,
)

data

# generate Weibull distribution
ref_data = {}
for n in data["Name"]:
    ref_data[n] = {}
    ref_data[n]["wind_speed"] = [0 + 0.01 * n for n in range(3100)]
    ref_data[n]["power_curve"] = []
    ref_data[n][n] = []
    for v in ref_data[n]["wind_speed"]:
        ref_data[n]["power_curve"].append(opt.ref_power_curve(v=v))
        ref_data[n][n].append(
            opt.weibull_probability_distribution(
                v=v,
                k=data[data["Name"] == n][("k", "mean")].iloc[0],
                c=data[data["Name"] == n][("c", "mean")].iloc[0],
            )
        )
    ref_data[n] = pd.DataFrame(ref_data[n])

ref_data = pd.concat(ref_data.values(), axis=1).T.drop_duplicates().T

ref_data.head()

ax = ref_data.plot(
    x="wind_speed",
    y="power_curve",
    ylabel="Power [MW]",
    linewidth=3,
    color=sns.color_palette("flare", 256)[127],
    figsize=(12, 6),
    legend=False,
)
ax.set_xlabel("Wind speed [m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")
ax.set_ylabel("Power [MW]")
sns.despine()

ax = ref_data.drop(columns=["power_curve"]).plot(
    x="wind_speed",
    cmap="flare_r",
    figsize=(12, 6),
    linestyle="dashed",
    linewidth=2,
    alpha=0.85,
)
ax.set_xlabel("Wind speed [m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")
ax.set_ylabel(
    "Weibull probability distribution function "
    "[s m\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]"
)
sns.despine()

# ## Annual energy production [MWh]

# max wind farm capacity
data["capacity"] = [1300, 824, 500]

# number of 15 MW turbines, rounded down to the nearest integer
data["n_turbines"] = (data["capacity"] / 15).astype(int)

aep = []
for n in data["Name"]:
    aepwt = opt.annual_energy_production(
        n_turbines=data[data["Name"] == n]["n_turbines"].iloc[0],
        k=data[data["Name"] == n][("k", "mean")].iloc[0],
        c=data[data["Name"] == n][("c", "mean")].iloc[0],
    )
    aep.append(aepwt)

aep = pd.DataFrame(aep)
aep.columns = ["AEP", "integral", "abserr"]

data = pd.concat([data, aep], axis=1)

# ## Annual hydrogen production [kg]

data["AHP"] = opt.annual_hydrogen_production(aep=data["AEP"])

# ## AHP as a proportion of the total working mass

data["AHP_frac"] = data["AHP"] / caverns[["working_mass"]].sum().iloc[0]

data

# ## Number of caverns required based on cumulative working mass and AHP

working_mass_cumsum_1 = (
    caverns.sort_values("working_mass", ascending=False)
    .reset_index()["working_mass"]
    .cumsum()
)

working_mass_cumsum_2 = (
    caverns.sort_values("working_mass").reset_index()["working_mass"].cumsum()
)

for x in range(len(data)):
    print(data["Name"].iloc[x])
    print("Working mass [kg]:", data["AHP"].iloc[x])
    print(
        "Number of caverns required:",
        working_mass_cumsum_1.loc[working_mass_cumsum_1 >= data["AHP"].iloc[x]]
        .head(1)
        .index[0]
        + 1,
        "-",
        working_mass_cumsum_2.loc[working_mass_cumsum_2 >= data["AHP"].iloc[x]]
        .head(1)
        .index[0]
        + 1,
    )
    print("-" * 50)

# ## Calculate distance between caverns and the wind farms and injection point [km]

wind_farms["Name_"] = wind_farms["Name"].str.split(expand=True)[0]
wind_farms = wind_farms.dissolve("Name_").reset_index()

# Dublin Port coordinates (Dinh et al. - injection point)
injection_point = gpd.GeoSeries(
    [Point(-(6 + 12 / 60), 53 + 21 / 60)], crs=4326
)
injection_point = injection_point.to_crs(rd.CRS)

distance_ip = []
for j in range(len(caverns)):
    distance_ip.append(
        injection_point.distance(
            caverns.iloc[[j]]["geometry"], align=False
        ).values[0]
        / 1000
    )
caverns["distance_ip"] = distance_ip

distance_wf = {}
for i in range(len(wind_farms)):
    distance_wf[wind_farms["Name_"][i]] = []
    for j in range(len(caverns)):
        distance_wf[wind_farms["Name_"][i]].append(
            (
                wind_farms.iloc[[i]]
                .distance(caverns.iloc[[j]]["geometry"], align=False)
                .values[0]
                / 1000
                + caverns.iloc[[j]]["distance_ip"].values[0]
            )
        )
    caverns[f"distance{wind_farms['Name_'][i]}"] = distance_wf[
        wind_farms["Name_"][i]
    ]

caverns = caverns.rename(columns={"distanceNorth": "distanceNISA"})

# ## CAPEX for pipeline [€ km⁻¹]

# 1,000 MW electrolyser
capex = opt.capex_pipeline(e_cap=1000)

capex

# ## LCOT for pipeline [€ kg⁻¹]

for wf in ["Codling", "Dublin", "NISA"]:
    caverns[f"LCOT_{wf}"] = opt.lcot_pipeline(
        capex=capex,
        transmission_distance=caverns[f"distance{wf}"],
        ahp=data[data["Name"].str.contains(wf)]["AHP"].values[0],
    )

caverns[
    [
        "cavern_depth",
        "working_mass",
        "capacity",
        "distance_ip",
        "distanceCodling",
        "distanceDublin",
        "distanceNISA",
        "LCOT_Codling",
        "LCOT_Dublin",
        "LCOT_NISA",
    ]
].describe()

# ## Maps

shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)


def plot_map_facet(cavern_df, classes, fontsize=11.5):
    """Helper function for plotting LCOT facet maps"""
    fig = plt.figure(figsize=(11, 11.75))
    xmin_, ymin_, xmax_, ymax_ = cavern_df.total_bounds
    colours = [
        int(n * 255 / (len(classes) - 2)) for n in range(len(classes) - 1)
    ]
    legend_handles = []

    for a, wf in enumerate(["Codling", "Dublin", "NISA"]):
        ax = fig.add_subplot(2, 2, a + 1, projection=ccrs.epsg(rd.CRS))
        gpd.GeoDataFrame(cavern_df, geometry=cavern_df.centroid).plot(
            ax=ax,
            scheme="UserDefined",
            classification_kwds={"bins": classes[1:]},
            column=f"LCOT_{wf}",
            zorder=2,
            marker=".",
            cmap="flare",
            markersize=8,
        )
        shape.plot(
            ax=ax, color="white", alpha=0.5, edgecolor="slategrey", zorder=1
        )
        cx.add_basemap(
            ax,
            crs=rd.CRS,
            source=cx.providers.CartoDB.VoyagerNoLabels,
            attribution=False,
        )
        if a in (0, 2):
            draw_labels = {"bottom": "x", "left": "y"}
        else:
            draw_labels = {"bottom": "x"}
        ax.gridlines(
            draw_labels=draw_labels,
            color="none",
            xlabel_style={"fontsize": fontsize},
            ylabel_style={"fontsize": fontsize},
            xformatter=LongitudeFormatter(number_format=".2f"),
            yformatter=LatitudeFormatter(number_format=".2f"),
        )
        if a == 2:
            ax.add_artist(
                ScaleBar(
                    1,
                    box_alpha=0,
                    location="lower right",
                    color="darkslategrey",
                    font_properties={"size": fontsize},
                )
            )
        plt.xlim(xmin_ - 1000, xmax_ + 1000)
        plt.ylim(ymin_ - 1000, ymax_ + 1000)
        if wf == "Codling":
            ax.set_title(f"{wf} Wind Park")
        elif wf == "Dublin":
            ax.set_title(f"{wf} Array")
        else:
            ax.set_title("North Irish Sea Array")

    for n, c in enumerate(colours):
        if n == 0:
            l = f"< {classes[1:][n]}"
        elif n == len(colours) - 1:
            l = f"≥ {classes[1:][-2]}"
        else:
            l = f"{classes[1:][n - 1]} - {classes[1:][n]}"
        legend_handles.append(
            mpatches.Patch(
                facecolor=sns.color_palette("flare", 256)[c], label=l
            )
        )
        plt.legend(
            loc="lower right",
            bbox_to_anchor=(1, 0.075),
            handles=legend_handles,
            title="Pipeline LCOT [€ kg⁻¹]",
            fontsize=fontsize,
            title_fontsize=fontsize,
        )
    plt.tight_layout()
    plt.show()


plot_map_facet(caverns, [0.04 * n for n in range(6)])
