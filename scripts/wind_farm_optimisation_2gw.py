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
import seaborn as sns
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import capacity as cap
from h2ss import compare
from h2ss import data as rd
from h2ss import functions as fns
from h2ss import optimisation as opt

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

opt.hydrogen_pipeline_density(pressure=65e5, temperature=10)

opt.hydrogen_pipeline_density(pressure=100e5, temperature=10)

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
        "data", "wind-farms", "marine-area-consent-wind.zip"
    )
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
    data_path=os.path.join("data", "subsea-cables", "KIS-ORCA.gpkg"),
    dat_extent=extent,
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
weibull_wf_df = fns.read_weibull_data(
    data_path_weibull=os.path.join(
        "data", "weibull-parameters-wind-speeds", "Weibull_150m_params_ITM.zip"
    ),
    data_path_wind_farms=os.path.join(
        "data", "wind-farms", "marine-area-consent-wind.zip"
    ),
)

weibull_powercurve = opt.weibull_distribution(weibull_wf_data=weibull_wf_df)

# ## Number of reference wind turbines

# 2 GW of offshore wind for green hydrogen production by 2030
# in addition to 5 GW offshore wind target in CLimate Action Plan 2023
# pg. 134
# assume this 2 GW is distributed evenly to the total capacity
print(
    f"{2 / 7 * 100:.2f}% of total offshore wind farm capacity dedicated "
    f"to H\N{SUBSCRIPT TWO} production"
)

# max wind farm capacity
weibull_wf_df["cap"] = [int(x * 2 / 7) for x in [1300, 824, 500]]

# number of 15 MW turbines, rounded down to the nearest integer
weibull_wf_df["n_turbines"] = opt.number_of_turbines(
    owf_cap=weibull_wf_df["cap"]
)

# ## Annual energy production [MWh]

weibull_wf_df = opt.annual_energy_production(weibull_wf_data=weibull_wf_df)

# ## Annual hydrogen production [kg]

weibull_wf_df["AHP"] = opt.annual_hydrogen_production(aep=weibull_wf_df["AEP"])

# ## AHP as a proportion of the total working mass

weibull_wf_df["AHP_frac"] = (
    weibull_wf_df["AHP"] / caverns[["working_mass"]].sum().iloc[0]
)

# ## AHP converted to storage demand [GWh]

weibull_wf_df["demand"] = cap.energy_storage_capacity(
    m_working=weibull_wf_df["AHP"]
)

# ## Number of caverns required based on cumulative working mass and AHP

compare.calculate_number_of_caverns(
    cavern_df=caverns, weibull_wf_data=weibull_wf_df
)

# ## Transmission distance [km]

caverns, injection_point = opt.transmission_distance(
    cavern_df=caverns, wf_data=wind_farms
)

# ## Electrolyser capacity [MW]

weibull_wf_df["E_cap"] = opt.electrolyser_capacity(
    n_turbines=weibull_wf_df["n_turbines"]
)

# ## CAPEX for pipeline [€ km⁻¹]

weibull_wf_df["CAPEX"] = opt.capex_pipeline(e_cap=weibull_wf_df["E_cap"])

weibull_wf_df

# totals
weibull_wf_df[
    ["cap", "n_turbines", "AEP", "AHP", "AHP_frac", "demand", "E_cap", "CAPEX"]
].sum()

compare.electricity_demand_ie(data=weibull_wf_df["demand"])

compare.hydrogen_demand_ie(data=weibull_wf_df["demand"])

# ## LCOT for pipeline [€ kg⁻¹]

caverns = opt.lcot_pipeline(weibull_wf_data=weibull_wf_df, cavern_df=caverns)

caverns[
    [
        "cavern_depth",
        "working_mass",
        "capacity",
        "distance_ip",
    ]
    + list(caverns.filter(like="dist_"))
    + list(caverns.filter(like="LCOT_"))
].describe()

caverns[list(caverns.filter(like="LCOT_"))].describe().mean(axis=1)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
sns.boxplot(
    caverns.filter(like="dist_")
    .set_axis(list(wind_farms["name"]), axis=1)
    .melt(),
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
    caverns.filter(like="LCOT_")
    .set_axis(list(wind_farms["name"]), axis=1)
    .melt(),
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
#     os.path.join("graphics", "fig_box_transmission_ntg_2gw.jpg"),
#     format="jpg",
#     dpi=600,
# )
plt.show()

# ## Maps

shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)


def plot_map_facet(
    cavern_df,
    classes,
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
            label = f"{classes[n1 - 1]:.2f}–<{classes[n1]:.2f}"
        legend_handles.append(
            mpatches.Patch(
                facecolor=sns.color_palette("flare", 256)[c], label=label
            )
        )

    for a, wf1 in enumerate(list(wind_farms["name"])):
        ax1 = fig1.add_subplot(2, 2, a + 1, projection=ccrs.epsg(rd.CRS))
        gpd.GeoDataFrame(cavern_df, geometry=cavern_df.centroid).plot(
            ax=ax1,
            scheme="UserDefined",
            classification_kwds={"bins": classes},
            column=f"LCOT_{wf1.replace(' ', '_')}",
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
                ylabel_style={"fontsize": fontsize, "rotation": 89.9},
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
        ax1.set_title(list(wind_farms["name"])[a])

    plt.tight_layout()
    # plt.savefig(
    #     os.path.join("graphics", "fig_map_transmission_ntg_2gw.jpg"),
    #     format="jpg", dpi=600
    # )
    plt.show()


plot_map_facet(
    caverns,
    list(mc.Quantiles(caverns[list(caverns.filter(like="LCOT_"))], k=6).bins),
)
