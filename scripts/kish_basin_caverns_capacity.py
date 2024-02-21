#!/usr/bin/env python
# coding: utf-8

# # Cavern storage capacity

import os

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import ticker
from matplotlib.lines import Line2D
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import capacity as cap
from h2ss import data as rd
from h2ss import functions as fns

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# ## Halite data

ds, extent = rd.kish_basin_data_depth_adjusted(
    dat_path=os.path.join("data", "kish-basin"),
    bathymetry_path=os.path.join("data", "bathymetry"),
)

xmin, ymin, xmax, ymax = extent.total_bounds

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

# height = 85 m, 500 m <= depth <= 2,000 m, diameter = 80 m,
# separation = 320 m
zones, zds = fns.zones_of_interest(
    dat_xr=ds, constraints={"height": 85, "min_depth": 500, "max_depth": 2000}
)

# ## Generate caverns

caverns, _ = fns.generate_caverns_with_constraints(
    zones_gdf=zones,
    zones_ds=zds,
    dat_extent=extent,
    depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
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
    heights=[85, 155, 311],
    depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
)

# ## Capacity

# ### Volume

caverns["cavern_volume"] = cap.cavern_volume(height=caverns["cavern_height"])
caverns["cavern_volume"] = cap.corrected_cavern_volume(
    v_cavern=caverns["cavern_volume"]
)

# ### Mid-point temperature

caverns["t_mid_point"] = cap.temperature_cavern_mid_point(
    height=caverns["cavern_height"], depth_top=caverns["cavern_depth"]
)

# ### Operating pressure

(
    caverns["p_operating_min"],
    caverns["p_operating_max"],
) = cap.pressure_operating(thickness_overburden=caverns["TopDepthSeabed"])

# ### Hydrogen gas density

caverns["rho_min"], caverns["rho_max"] = cap.density_hydrogen_gas(
    p_operating_min=caverns["p_operating_min"],
    p_operating_max=caverns["p_operating_max"],
    t_mid_point=caverns["t_mid_point"],
)

# ### Working mass of hydrogen

(
    caverns["working_mass"],
    caverns["mass_operating_min"],
    caverns["mass_operating_max"],
) = cap.mass_hydrogen_working(
    rho_h2_min=caverns["rho_min"],
    rho_h2_max=caverns["rho_max"],
    v_cavern=caverns["cavern_volume"],
)

# ### Energy storage capacity in GWh

caverns["capacity"] = cap.energy_storage_capacity(
    m_working=caverns["working_mass"]
)

# ## Stats

# proportion of working gas to total gas
caverns["working_mass_pct"] = caverns["working_mass"] / (
    caverns["working_mass"] + caverns["mass_operating_min"]
)

caverns.drop(
    ["x", "y", "TopTWT", "BaseDepth", "TopDepth", "BaseDepthSeabed"], axis=1
).describe()

# cavern volumes
list(caverns["cavern_volume"].unique())

# totals
caverns[
    [
        "cavern_volume",
        "working_mass",
        "capacity",
        "mass_operating_min",
        "mass_operating_max",
    ]
].sum()

# compare to Ireland's electricity demand in 2050 (Deane, 2021)
print(
    "Energy capacity as a percentage of Ireland's electricity demand in 2050:",
    f"{(caverns['capacity'].sum() / 1000 / 122 * 100):.2f}–"
    + f"{(caverns['capacity'].sum() / 1000 / 84 * 100):.2f}%",
)

# total capacity at various depth/height combinations
s = caverns.groupby(["height", "depth"], sort=False)[["capacity"]].sum()
s["%"] = s["capacity"] / caverns[["capacity"]].sum().iloc[0] * 100
s

s.groupby("depth").sum()[["capacity"]]

s.groupby("height").sum()[["capacity"]]

# number of caverns
s = caverns.groupby(["height", "depth"], sort=False)[["capacity"]].count()
s["%"] = s["capacity"] / len(caverns) * 100
s

s.groupby("depth").sum()[["capacity"]]

s.groupby("height").sum()[["capacity"]]

fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
sns.histplot(
    caverns.rename(columns={"depth": "Cavern top depth [m]"}).sort_values(
        "TopDepthSeabed"
    ),
    x="capacity",
    hue="Cavern top depth [m]",
    palette="rocket_r",
    multiple="stack",
    alpha=1,
    ax=axes[0],
)
sns.histplot(
    caverns.rename(columns={"height": "Cavern height [m]"}).sort_values(
        "Thickness"
    ),
    x="capacity",
    hue="Cavern height [m]",
    palette="rocket_r",
    multiple="stack",
    alpha=1,
    ax=axes[1],
)
axes[0].set_xlabel("Hydrogen storage capacity [GWh]")
axes[1].set_xlabel("Hydrogen storage capacity [GWh]")
axes[0].grid(which="major", axis="y")
axes[1].grid(which="major", axis="y")
axes[0].set_ylabel("Number of caverns")
sns.despine()
plt.tight_layout()
plt.show()

# copy dataframe
caverns_pot_all = caverns.copy()

# ## Maps

# create exclusion buffer
buffer = pd.concat([wells_b, shipwrecks_b, shipping_b, cables_b]).dissolve()


def plot_map_alt(
    dat_xr, cavern_df, zones_gdf, classes, top_depth=True, fontsize=11.5
):
    """Helper function to plot caverns within the zones of interest"""
    plt.figure(figsize=(20, 11.5))
    axis1 = plt.axes(projection=ccrs.epsg(rd.CRS))
    legend_handles1 = []
    classes = sorted(classes)

    # halite boundary - use buffering to smooth the outline
    shape = rd.halite_shape(dat_xr=dat_xr).buffer(1000).buffer(-1000)
    shape.plot(
        ax=axis1,
        edgecolor="darkslategrey",
        color="none",
        linewidth=2,
        alpha=0.5,
        zorder=2,
    )
    legend_handles1.append(
        mpatches.Patch(
            facecolor="none",
            linewidth=2,
            edgecolor="darkslategrey",
            label="Kish Basin boundary",
            alpha=0.5,
        )
    )

    zones_gdf.plot(
        ax=axis1, zorder=1, linewidth=0, facecolor="white", alpha=0.45
    )
    zones_gdf.plot(
        ax=axis1,
        zorder=2,
        edgecolor="slategrey",
        linestyle="dotted",
        linewidth=1.25,
        facecolor="none",
    )
    legend_handles1.append(
        mpatches.Patch(
            facecolor="none",
            linestyle="dotted",
            edgecolor="slategrey",
            label="Area of interest",
            linewidth=1.25,
        )
    )

    pd.concat([buffer, wind_farms]).dissolve().clip(shape).plot(
        ax=axis1,
        facecolor="none",
        linewidth=0.65,
        edgecolor="slategrey",
        zorder=2,
        alpha=0.5,
        hatch="//",
    )
    legend_handles1.append(
        mpatches.Patch(
            facecolor="none",
            hatch="//",
            edgecolor="slategrey",
            label="Exclusion",
            alpha=0.65,
            linewidth=0.5,
        )
    )

    legend_handles1.append(
        mpatches.Patch(
            label="Hydrogen storage \ncapacity [GWh]", visible=False
        )
    )

    colours = [int(n * 255 / (len(classes) - 1)) for n in range(len(classes))]
    for n, y in enumerate(colours):
        if n == 0:
            c = cavern_df[cavern_df["capacity"] < classes[1]]
            label1 = f"< {classes[1]}"
        elif n == len(classes) - 1:
            c = cavern_df[cavern_df["capacity"] >= classes[n]]
            label1 = f"≥ {classes[n]}"
        else:
            c = cavern_df[
                (cavern_df["capacity"] >= classes[n])
                & (cavern_df["capacity"] < classes[n + 1])
            ]
            label1 = f"{classes[n]}–{classes[n + 1]}"
        if top_depth:
            for df, markersize in zip(
                [
                    c[c["depth"] == "500 - 1,000"],
                    c[c["depth"] == "1,000 - 1,500"],
                    c[c["depth"] == "1,500 - 2,000"],
                ],
                [20, 50, 20],
            ):
                if len(df) > 0:
                    df.centroid.plot(
                        ax=axis1,
                        zorder=3,
                        linewidth=0,
                        marker=".",
                        markersize=markersize,
                        color=sns.color_palette("flare", 256)[y],
                    )
        else:
            gpd.GeoDataFrame(cavern_df, geometry=cavern_df.centroid).plot(
                ax=axis1,
                scheme="UserDefined",
                classification_kwds={"bins": classes[1:]},
                column="capacity",
                zorder=3,
                marker=".",
                cmap="flare",
                markersize=20,
            )
        legend_handles1.append(
            mpatches.Patch(
                facecolor=sns.color_palette("flare", 256)[y], label=label1
            )
        )

    if top_depth:
        legend_handles1.append(
            mpatches.Patch(label="Cavern top depth [m]", visible=False)
        )
        for markersize, label1 in zip(
            [6, 3], ["1,000–1,500", "500–1,000 or \n1,500–2,000"]
        ):
            legend_handles1.append(
                Line2D(
                    [0],
                    [0],
                    marker=".",
                    linewidth=0,
                    label=label1,
                    color="darkslategrey",
                    markersize=markersize,
                )
            )

    plt.xlim(shape.bounds["minx"][0] - 1000, shape.bounds["maxx"][0] + 1000)
    plt.ylim(shape.bounds["miny"][0] - 1000, shape.bounds["maxy"][0] + 1000)

    cx.add_basemap(
        axis1, crs=rd.CRS, source=cx.providers.CartoDB.VoyagerNoLabels, zoom=12
    )
    axis1.gridlines(
        draw_labels={"bottom": "x", "left": "y"},
        alpha=0.25,
        color="darkslategrey",
        xlabel_style={"fontsize": fontsize},
        ylabel_style={"fontsize": fontsize},
    )
    axis1.add_artist(
        ScaleBar(
            1,
            box_alpha=0,
            location="lower right",
            color="darkslategrey",
            width_fraction=0.0075,
            font_properties={"size": fontsize},
        )
    )
    plt.legend(
        loc="lower right",
        bbox_to_anchor=(1, 0.05),
        handles=legend_handles1,
        fontsize=fontsize,
    )

    plt.tight_layout()
    plt.show()


plot_map_alt(ds, caverns, zones, [40 * n for n in range(5)])

# ## Restrict cavern height to 155 m, depth to 1,000-1,500 m

# height = 155 m, 1,000 m <= depth <= 1,500 m, diameter = 80 m,
# separation = 320 m
zones, zds = fns.zones_of_interest(
    dat_xr=ds,
    constraints={"height": 155, "min_depth": 1000, "max_depth": 1500},
)

caverns, _ = fns.generate_caverns_with_constraints(
    zones_gdf=zones,
    zones_ds=zds,
    dat_extent=extent,
    depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
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

# calculate volumes and capacities
caverns = cap.calculate_capacity_dataframe(cavern_df=caverns)

# proportion of working gas to total gas
caverns["working_mass_pct"] = caverns["working_mass"] / (
    caverns["working_mass"] + caverns["mass_operating_min"]
)

caverns.drop(
    [
        "x",
        "y",
        "TopTWT",
        "BaseDepth",
        "TopDepth",
        "BaseDepthSeabed",
        "height",
        "cavern_height",
    ],
    axis=1,
).describe()

# totals
caverns[
    [
        "cavern_volume",
        "working_mass",
        "capacity",
        "mass_operating_min",
        "mass_operating_max",
    ]
].sum()

# compare to Ireland's electricity demand in 2050 (Deane, 2021)
print(
    "Energy capacity as a percentage of Ireland's electricity demand in 2050:",
    f"{(caverns['capacity'].sum() / 1000 / 122 * 100):.2f}–"
    + f"{(caverns['capacity'].sum() / 1000 / 84 * 100):.2f}%",
)

plot_map_alt(ds, caverns, zones, [80 + n * 5 for n in range(6)], False)

# ## Distribution


def cavern_boxplot(cavern_df):
    """Helper function for creating boxplots"""
    _, axes1 = plt.subplots(1, 4, figsize=(11, 4.5))
    sns.boxplot(
        cavern_df,
        y="cavern_depth",
        color=sns.color_palette("rocket", 1)[0],
        width=0.2,
        ax=axes1[0],
        legend=False,
        linecolor="black",
        linewidth=1.1,
        gap=0.1,
        flierprops={"markeredgecolor": "grey", "alpha": 0.5},
    )
    axes1[0].set_ylabel("Top depth [m]")
    axes1[0].get_yaxis().set_major_formatter(
        ticker.FuncFormatter(lambda x, p: format(int(x), ","))
    )
    sns.boxplot(
        (cavern_df[["p_operating_min", "p_operating_max"]] / 1e6)
        .rename(
            columns={
                "p_operating_min": "min",
                "p_operating_max": "max",
            }
        )
        .melt(),
        y="value",
        hue="variable",
        palette="rocket_r",
        width=0.4,
        ax=axes1[1],
        linecolor="black",
        linewidth=1.1,
        gap=0.1,
        flierprops={"markeredgecolor": "grey", "alpha": 0.5},
    )
    axes1[1].set_ylabel("Operating pressure [MPa]")
    axes1[1].legend()
    sns.boxplot(
        (cavern_df[["mass_operating_min", "working_mass"]] / 1e6)
        .rename(
            columns={
                "mass_operating_min": "cushion",
                "working_mass": "working",
            }
        )
        .melt(),
        y="value",
        hue="variable",
        palette="rocket_r",
        width=0.4,
        ax=axes1[2],
        linecolor="black",
        linewidth=1.1,
        gap=0.1,
        flierprops={"markeredgecolor": "grey", "alpha": 0.5},
    )
    axes1[2].set_ylabel("Gas mass [kt]")
    axes1[2].legend()
    sns.boxplot(
        cavern_df,
        y="capacity",
        color=sns.color_palette("rocket", 1)[0],
        width=0.2,
        ax=axes1[3],
        legend=False,
        linecolor="black",
        linewidth=1.1,
        gap=0.1,
        flierprops={"markeredgecolor": "grey", "alpha": 0.5},
    )
    axes1[3].set_ylabel("Energy storage capacity [GWh]")
    sns.despine(bottom=True)
    plt.tight_layout()
    plt.show()


cavern_boxplot(caverns_pot_all)

cavern_boxplot(caverns)

fig, axes = plt.subplots(3, 2, figsize=(12, 8))
for variable, label, axis in zip(
    [
        "cavern_depth",
        "working_mass",
        "p_operating_min",
        "mass_operating_min",
        "p_operating_max",
        "capacity",
    ],
    [
        "Top depth [m]",
        "Working gas mass [kt]",
        "Minimum operating pressure [MPa]",
        "Cushion gas mass [kt]",
        "Maximum operating pressure [MPa]",
        "Energy storage capacity [GWh]",
    ],
    axes.flat,
):
    if variable in ["cavern_depth", "capacity"]:
        d1 = caverns_pot_all[[variable]]
        d2 = caverns[[variable]]
    else:
        d1 = caverns_pot_all[[variable]] / 1e6
        d2 = caverns[[variable]] / 1e6
    sns.boxplot(
        pd.concat(
            [d1.set_axis(["all"], axis=1), d2.set_axis(["optimal"], axis=1)]
        )
        .melt()
        .dropna(),
        x="value",
        hue="variable",
        palette="flare_r",
        linecolor="black",
        linewidth=1.1,
        gap=0.3,
        width=0.55,
        flierprops={"markeredgecolor": "grey", "alpha": 0.5},
        ax=axis,
        legend=False,
    )
    axis.set_xlabel(label)
    if variable == "cavern_depth":
        axis.get_xaxis().set_major_formatter(
            ticker.FuncFormatter(lambda x, p: format(int(x), ","))
        )

legend_handles = [
    mpatches.Patch(
        facecolor=sns.color_palette("flare_r", 2)[0],
        label="85–311 m tall caverns at 500–2,000 m depth",
        edgecolor="black",
    ),
    mpatches.Patch(
        facecolor=sns.color_palette("flare_r", 2)[1],
        label="155 m tall caverns at 1,000–1,500 m depth",
        edgecolor="black",
    ),
]
plt.legend(loc="lower right", handles=legend_handles, fontsize=11.5)
sns.despine()
plt.tight_layout()
plt.show()
