#!/usr/bin/env python
# coding: utf-8

# # Wind farm optimisation - 2 GW of dedicated offshore wind for hydrogen production

import os

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import mapclassify as mc
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib import patheffects
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

zones, zds = fns.zones_of_interest(
    dat_xr=ds,
    constraints={"net_height": 120, "min_depth": 500, "max_depth": 2000},
)

# ## Generate caverns

caverns = fns.generate_caverns_hexagonal_grid(
    zones_df=zones,
    dat_extent=extent,
)

caverns = fns.cavern_dataframe(
    dat_zone=zds,
    cavern_df=caverns,
    depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
)

# label caverns by depth and heights
caverns = fns.label_caverns(
    cavern_df=caverns,
    heights=[120],
    depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
)

caverns, _ = fns.generate_caverns_with_constraints(
    cavern_df=caverns,
    exclusions={
        "wells": wells_b,
        "wind_farms": wind_farms,
        "shipwrecks": shipwrecks_b,
        "shipping": shipping_b,
        "cables": cables_b,
        "edge": edge_buffer,
    },
)

# ## Capacity

caverns["cavern_total_volume"] = cap.cavern_volume(
    height=caverns["cavern_height"]
)
caverns["cavern_volume"] = cap.corrected_cavern_volume(
    v_cavern=caverns["cavern_total_volume"]
)

caverns["t_mid_point"] = cap.temperature_cavern_mid_point(
    height=caverns["cavern_height"], depth_top=caverns["cavern_depth"]
)

(
    caverns["p_operating_min"],
    caverns["p_operating_max"],
) = cap.pressure_operating(
    thickness_overburden=caverns["TopDepthSeabed"],
    depth_water=-caverns["Bathymetry"],
)

caverns["rho_min"], caverns["rho_max"] = cap.density_hydrogen_gas(
    p_operating_min=caverns["p_operating_min"],
    p_operating_max=caverns["p_operating_max"],
    t_mid_point=caverns["t_mid_point"],
)

(
    caverns["working_mass"],
    caverns["mass_operating_min"],
    caverns["mass_operating_max"],
) = cap.mass_hydrogen_working(
    rho_h2_min=caverns["rho_min"],
    rho_h2_max=caverns["rho_max"],
    v_cavern=caverns["cavern_volume"],
)

caverns["capacity"] = cap.energy_storage_capacity(
    m_working=caverns["working_mass"]
)

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

# generate Weibull distribution
ref_data = {}
for n in data["Name"]:
    ref_data[n] = {}
    ref_data[n]["wind_speed"] = [0 + 0.01 * n for n in range(3000)]
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

# ## Annual energy production [MWh]

# 2 GW of offshore wind for green hydrogen production by 2030
# in addition to 5 GW offshore wind target in CLimate Action Plan 2023
# pg. 134
# assume this 2 GW is distributed evenly to the total capacity
print(
    f"{2000 / 7000 * 100:.2f}% of total offshore wind farm capacity dedicated",
    f"to for H\N{SUBSCRIPT TWO} production",
)

# max wind farm capacity
data["cap"] = [int(x * 2000 / 7000) for x in [1300, 824, 500]]

# number of 15 MW turbines, rounded down to the nearest integer
data["n_turbines"] = opt.number_of_turbines(owf_cap=data["cap"])

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

# ## AHP converted to storage demand [GWh]

data["demand"] = cap.energy_storage_capacity(m_working=data["AHP"])

# ## Number of caverns required based on cumulative working mass and AHP

working_mass_cumsum_1 = (
    caverns.sort_values("working_mass", ascending=False)
    .reset_index()[["working_mass", "capacity"]]
    .cumsum()
)

working_mass_cumsum_2 = (
    caverns.sort_values("working_mass")
    .reset_index()[["working_mass", "capacity"]]
    .cumsum()
)

caverns_low = []
caverns_high = []
cap_max = []
for x in range(len(data)):
    print(data["Name"].iloc[x])
    print(f"Working mass [kg]: {(data['AHP'].iloc[x]):.6E}")
    caverns_low.append(
        working_mass_cumsum_1.loc[
            working_mass_cumsum_1["working_mass"] >= data["AHP"].iloc[x]
        ]
        .head(1)
        .index[0]
        + 1
    )
    caverns_high.append(
        working_mass_cumsum_2.loc[
            working_mass_cumsum_2["working_mass"] >= data["AHP"].iloc[x]
        ]
        .head(1)
        .index[0]
        + 1
    )
    print(f"Number of caverns required: {caverns_low[x]}–{caverns_high[x]}")
    cap_max.append(
        max(
            working_mass_cumsum_1.loc[
                working_mass_cumsum_1["working_mass"] >= data["AHP"].iloc[x]
            ]
            .head(1)["capacity"]
            .values[0],
            working_mass_cumsum_2.loc[
                working_mass_cumsum_2["working_mass"] >= data["AHP"].iloc[x]
            ]
            .head(1)["capacity"]
            .values[0],
        )
    )
    print(f"Capacity (approx.) [GWh]: {(cap_max[x]):.2f}")
    print("-" * 50)

# total number of caverns
print(
    f"Total number of caverns required: {sum(caverns_low)}–{sum(caverns_high)}"
)

# number of caverns as a percentage of the total available caverns
print(
    "Number of caverns required as a percentage of all available caverns:",
    f"{(sum(caverns_low) / len(caverns) * 100):.2f}–"
    + f"{(sum(caverns_high) / len(caverns) * 100):.2f}%",
)

# total capacity
print(f"Total maximum cavern capacity (approx.): {sum(cap_max):.2f} GWh")

# ## Transmission distance [km]

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
    caverns[f"dist_{wind_farms['Name_'][i]}"] = distance_wf[
        wind_farms["Name_"][i]
    ]

# ## CAPEX for pipeline [€ km⁻¹]

# calculate electrolyser capacity
data["E_cap"] = opt.electrolyser_capacity(n_turbines=data["n_turbines"])

data["CAPEX"] = opt.capex_pipeline(e_cap=data["E_cap"])

data

# totals
data[
    ["cap", "n_turbines", "AEP", "AHP", "AHP_frac", "demand", "E_cap", "CAPEX"]
].sum()

# compare to Ireland's electricity demand in 2050 (Deane, 2021)
print(
    "Storage demand as a percentage of Ireland's electricity demand in 2050:",
    f"{(sum(data['demand']) / 1000 / 122 * 100):.2f}–"
    + f"{(sum(data['demand']) / 1000 / 84 * 100):.2f}%",
)

# ## LCOT for pipeline [€ kg⁻¹]

for wf in list(wind_farms["Name_"]):
    caverns[f"LCOT_{wf}"] = opt.lcot_pipeline(
        capex=data[data["Name"].str.contains(wf)]["CAPEX"].values[0],
        transmission_distance=caverns[f"dist_{wf}"],
        ahp=data[data["Name"].str.contains(wf)]["AHP"].values[0],
    )

caverns[
    [
        "cavern_depth",
        "working_mass",
        "capacity",
        "distance_ip",
        "dist_Codling",
        "dist_Dublin",
        "dist_North",
        "LCOT_Codling",
        "LCOT_Dublin",
        "LCOT_North",
    ]
].describe()

caverns[
    [
        "LCOT_Codling",
        "LCOT_Dublin",
        "LCOT_North",
    ]
].describe().mean(axis=1)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
sns.boxplot(
    caverns.filter(like="dist_").set_axis(list(data["Name"]), axis=1).melt(),
    y="value",
    hue="variable",
    palette=sns.color_palette(["tab:red", "tab:gray", "tab:blue"]),
    width=0.35,
    ax=axes[0],
    legend=False,
    linecolor="black",
    linewidth=1.1,
    gap=0.15,
)
axes[0].set_ylabel("Transmission distance [km]")
axes[0].tick_params(axis="x", bottom=False)
sns.boxplot(
    caverns.filter(like="LCOT_").set_axis(list(data["Name"]), axis=1).melt(),
    y="value",
    hue="variable",
    palette=sns.color_palette(["tab:red", "tab:gray", "tab:blue"]),
    width=0.35,
    ax=axes[1],
    linecolor="black",
    linewidth=1.1,
    gap=0.15,
)
axes[1].set_ylabel("Pipeline LCOT [€ kg⁻¹]")
axes[1].legend(loc="lower right")
axes[1].tick_params(axis="x", bottom=False)
axes[1].yaxis.grid(True, linewidth=0.25)
axes[0].yaxis.grid(True, linewidth=0.25)
sns.despine(bottom=True)
plt.tight_layout()
# plt.savefig(
#     os.path.join("graphics", "fig_box_transmission_ntg.jpg"),
#     format="jpg",
#     dpi=600,
# )
plt.show()

# ## Maps

shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)


def plot_map_facet(
    cavern_df,
    classes,
    filename="fig_map_transmission_ntg.jpg",
    fontsize=11.5,
):
    """Helper function for plotting LCOT facet maps"""
    fig1 = plt.figure(figsize=(11, 11.75))
    xmin_, ymin_, xmax_, ymax_ = cavern_df.total_bounds
    colours = [int(n * 255 / (len(classes) - 1)) for n in range(len(classes))]
    legend_handles = []
    classes = sorted(classes)

    for n1, c in enumerate(colours):
        if n1 == 0:
            label = f"< {classes[n1]:.2f}"
        elif n1 == len(colours) - 1:
            label = f"≥ {classes[-2]:.2f}"
        else:
            label = f"{classes[n1 - 1]:.2f}–{classes[n1]:.2f}"
        legend_handles.append(
            mpatches.Patch(
                facecolor=sns.color_palette("flare", 256)[c], label=label
            )
        )

    for a, wf1 in enumerate(list(wind_farms["Name_"])):
        ax1 = fig1.add_subplot(2, 2, a + 1, projection=ccrs.epsg(rd.CRS))
        gpd.GeoDataFrame(cavern_df, geometry=cavern_df.centroid).plot(
            ax=ax1,
            scheme="UserDefined",
            classification_kwds={"bins": classes},
            column=f"LCOT_{wf1}",
            zorder=2,
            marker=".",
            cmap="flare",
            markersize=20,
        )
        shape.plot(
            ax=ax1, color="white", alpha=0.5, edgecolor="slategrey", zorder=1
        )
        cx.add_basemap(
            ax1,
            crs=rd.CRS,
            source=cx.providers.CartoDB.VoyagerNoLabels,
            attribution=False,
        )
        ax1.gridlines(
            draw_labels={"bottom": "x"},
            color="lightslategrey",
            alpha=0.25,
            xlabel_style={"fontsize": fontsize},
            xformatter=LongitudeFormatter(auto_hide=False, dms=True),
        )
        if not a == 1:
            ax1.gridlines(
                draw_labels={"left": "y"},
                color="lightslategrey",
                alpha=0.25,
                ylabel_style={"fontsize": fontsize, "rotation": 90},
                yformatter=LatitudeFormatter(auto_hide=False, dms=True),
            )
        if a == 2:
            ax1.add_artist(
                ScaleBar(
                    1,
                    box_alpha=0,
                    location="lower right",
                    color="darkslategrey",
                    font_properties={"size": fontsize},
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
        plt.xlim(xmin_ - 1000, xmax_ + 1000)
        plt.ylim(ymin_ - 1000, ymax_ + 1000)
        ax1.set_title(list(data["Name"])[a])

    plt.tight_layout()
    # plt.savefig(
    #     os.path.join("graphics", filename),
    #     format="jpg",
    #     dpi=600,
    # )
    plt.show()


classes = mc.Quantiles(
    caverns[
        [
            "LCOT_Codling",
            "LCOT_Dublin",
            "LCOT_North",
        ]
    ],
    k=6,
)

plot_map_facet(caverns, list(classes.bins))
