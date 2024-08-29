"""Utility functions for plotting."""

import os

import branca.colormap as cm
import folium
import geopandas as gpd
import numpy as np
import pandas as pd
import seaborn as sns
from branca.element import MacroElement
from folium.plugins import (Fullscreen, GroupedLayerControl, MeasureControl,
                            MousePosition, StripePattern)
from jinja2 import Template

from h2ss import compare
from h2ss import data as rd


class BindColormap(MacroElement):
    """Binds a colormap to a given layer on Folium.

    Parameters
    ----------
    colormap : branca.colormap.ColorMap
        The colormap to bind.

    Notes
    -----
    Source:
    https://nbviewer.org/gist/BibMartin/f153aa957ddc5fadc64929abdee9ff2e
    """

    def __init__(self, layer, colormap):
        super(BindColormap, self).__init__()
        self.layer = layer
        self.colormap = colormap
        self._template = Template(
            """
        {% macro script(this, kwargs) %}
            {{this.colormap.get_name()}}.svg[0][0].style.display = 'block';
            {{this._parent.get_name()}}.on('overlayadd', function (eventLayer) {
                if (eventLayer.layer == {{this.layer.get_name()}}) {
                    {{this.colormap.get_name()}}.svg[0][0].style.display = 'block';
                }});
            {{this._parent.get_name()}}.on('overlayremove', function (eventLayer) {
                if (eventLayer.layer == {{this.layer.get_name()}}) {
                    {{this.colormap.get_name()}}.svg[0][0].style.display = 'none';
                }});
        {% endmacro %}
        """
        )  # noqa


def plot_interactive_map():
    """Plot an interactive map of the results using Folium."""
    ds, extent, exclusions = compare.load_all_data(keep_orig=True)

    caverns, zones, weibull_df, injection_point = (
        compare.optimisation_function(
            ds=ds,
            extent=extent,
            exclusions=exclusions,
            cavern_diameter=85,
            cavern_height=120,
        )
    )

    shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)

    col_list = [
        "geometry",
        "cavern_depth",
        "capacity",
        "Bathymetry",
        "NetToGross",
        "ThicknessNet",
        "dist_Codling_Wind_Park",
        "dist_Dublin_Array",
        "dist_North_Irish_Sea_Array",
        "LCOT_Codling_Wind_Park",
        "LCOT_Dublin_Array",
        "LCOT_North_Irish_Sea_Array",
    ]
    caverns_plt = caverns[col_list].copy()
    caverns_plt["LCOT_mean"] = caverns_plt[
        [
            "LCOT_Codling_Wind_Park",
            "LCOT_Dublin_Array",
            "LCOT_North_Irish_Sea_Array",
        ]
    ].mean(axis=1)
    for f in ["capacity", "LCOT_mean"]:
        caverns_plt[f"{f}_str"] = caverns_plt[f].copy()
        caverns_plt[f"{f}_str"] = [
            f"{x:,.3f}" for x in caverns_plt[f"{f}_str"]
        ]
    for f in col_list:
        if f == "NetToGross":
            caverns_plt[f] = [f"{x * 100:.0f}" for x in caverns[f]]
        elif f == "geometry":
            caverns_plt[f] = caverns_plt[f].buffer(100)
        elif f == "Bathymetry":
            caverns_plt[f] = [f"{-x:.3f}" for x in caverns[f]]
        elif f not in ["capacity", "LCOT_mean"]:
            caverns_plt[f] = [f"{x:,.3f}" for x in caverns[f]]

    # exclusions["wind_farms"] = pd.concat(
    #     [
    #         exclusions["wind_farms"][["geometry"]],
    #         weibull_df.drop(columns=["integral", "abserr"])
    #     ],
    #     axis=1
    # )
    exclusions["wind_farms"] = exclusions["wind_farms"][
        ["name", "geometry", "cap"]
    ]
    exclusions["wind_farms"]["cap"] = [
        f"{x:,.0f}" for x in exclusions["wind_farms"]["cap"]
    ]

    exclusions["shipwrecks"] = exclusions["shipwrecks"][
        ["VESSELNAME", "VESSELTYPE", "geometry"]
    ]
    exclusions["shipwrecks"]["VESSELNAME"] = exclusions["shipwrecks"][
        "VESSELNAME"
    ].fillna("Unknown")

    exclusions["wells"] = exclusions["wells"][
        ["RIG_NAME", "OPERATOR", "geometry"]
    ]

    exclusions["cables"] = exclusions["cables"].dissolve()[["geometry"]]

    exclusions["shipping"] = exclusions["shipping"].dissolve()[["geometry"]]

    # create exclusion buffer
    buffer = (
        pd.concat(
            [
                exclusions["wells_b"],
                exclusions["shipwrecks_b"],
                exclusions["shipping_b"],
                exclusions["cables_b"],
            ]
        )
        .dissolve()
        .overlay(gpd.GeoDataFrame(geometry=shape))
    )

    minx, miny, maxx, maxy = shape.to_crs(4326).bounds.values[0]
    m = folium.Map(
        tiles="cartodbvoyager",
        control_scale=True,
        location=[np.mean((miny, maxy)), np.mean((minx, maxx))],
        zoom_start=10,
        min_lat=miny,
        max_lat=maxy,
        min_lon=minx,
        max_lon=maxx,
    )

    # shipwrecks
    folium.GeoJson(
        exclusions["shipwrecks"].to_crs(4326),
        name="Shipwrecks",
        tooltip=folium.GeoJsonTooltip(
            fields=["VESSELNAME", "VESSELTYPE"],
            aliases=["Shipwreck", "Vessel type"],
        ),
        marker=folium.Marker(icon=folium.Icon(icon="sailboat", prefix="fa")),
        style_function=lambda feature: {"markerColor": "gray"},
    ).add_to(m)

    # exploration wells
    folium.GeoJson(
        exclusions["wells"].to_crs(4326),
        name="Exploration wells",
        marker=folium.Marker(icon=folium.Icon(icon="oil-well", prefix="fa")),
        style_function=lambda feature: {"markerColor": "gray"},
        tooltip=folium.GeoJsonTooltip(
            fields=["RIG_NAME", "OPERATOR"],
            aliases=["Exploration well", "Operator"],
        ),
    ).add_to(m)

    # shipping routes
    folium.GeoJson(
        exclusions["shipping"].to_crs(4326),
        name="Frequent shipping routes",
        color="crimson",
        tooltip="Frequent shipping route",
    ).add_to(m)

    # subsea cables
    folium.GeoJson(
        exclusions["cables"].to_crs(4326),
        name="Subsea cables",
        color="royalblue",
        tooltip="Subsea cable",
    ).add_to(m)

    # wind farms
    folium.GeoJson(
        exclusions["wind_farms"].to_crs(4326),
        name="Wind farms",
        tooltip=folium.GeoJsonTooltip(
            fields=["name", "cap"], aliases=["Wind farm", "Capacity [MW]"]
        ),
        style_function=lambda feature: {
            "fillColor": "none",
            "color": "seagreen",
            "fillPattern": StripePattern(
                angle=135, color="seagreen", space_color="none"
            ),
            "weight": 1,
        },
    ).add_to(m)

    # caverns - capacity
    fg1 = folium.FeatureGroup(name="Cavern capacity")
    cm1 = cm.StepColormap(
        list(sns.color_palette("rocket_r")),
        vmin=caverns_plt["capacity"].min(),
        vmax=caverns_plt["capacity"].max(),
    )
    caverns_dict = caverns_plt.reset_index().set_index("index")["capacity"]
    cd1 = {key: cm1(caverns_dict[key]) for key in caverns_dict.keys()}
    folium.GeoJson(
        caverns_plt.to_crs(4326),
        style_function=lambda feature: {
            "fillColor": cd1[int(feature["id"])],
            "color": "darkslategrey",
            "weight": 0.5,
            "fillOpacity": 1,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=[
                "capacity_str",
                "cavern_depth",
                "ThicknessNet",
                "NetToGross",
                "Bathymetry",
            ],
            aliases=[
                "Cavern capacity [GWh]",
                "Cavern top depth [m]",
                "Net halite thickness [m]",
                "Net-to-gross ratio [%]",
                "Water depth [m]",
            ],
        ),
        smooth_factor=0,
    ).add_to(fg1)
    m.add_child(fg1)
    cm1.caption = "Cavern capacity [GWh]"

    # caverns - LCOT
    fg2 = folium.FeatureGroup(name="Pipeline LCOT")
    cm2 = cm.StepColormap(
        list(sns.color_palette("mako_r")),
        vmin=caverns_plt["LCOT_mean"].min(),
        vmax=caverns_plt["LCOT_mean"].max(),
    )
    caverns_dict = caverns_plt.reset_index().set_index("index")["LCOT_mean"]
    cd2 = {key: cm2(caverns_dict[key]) for key in caverns_dict.keys()}
    folium.GeoJson(
        caverns_plt.to_crs(4326),
        style_function=lambda feature: {
            "fillColor": cd2[int(feature["id"])],
            "color": "darkslategrey",
            "weight": 0.5,
            "fillOpacity": 1,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=[
                "LCOT_mean_str",
                "LCOT_Codling_Wind_Park",
                "LCOT_Dublin_Array",
                "LCOT_North_Irish_Sea_Array",
            ],
            aliases=[
                "Pipeline LCOT [€ kg⁻¹]<br>Mean",
                "Codling Wind Park",
                "Dublin Array",
                "North Irish Sea Array",
            ],
        ),
        smooth_factor=0,
    ).add_to(fg2)
    m.add_child(fg2)
    cm2.caption = "Mean pipeline LCOT [€ kg⁻¹]"

    # basin boundary
    folium.GeoJson(
        shape.to_crs(4326),
        name="Kish Basin boundary",
        style_function=lambda feature: {
            "fillColor": "none",
            "color": "darkslategrey",
            "weight": 1,
        },
        tooltip="Kish Basin boundary",
    ).add_to(m)

    # exclusion buffer
    folium.GeoJson(
        buffer.to_crs(4326),
        name="Exclusion buffer",
        style_function=lambda feature: {
            "fillPattern": StripePattern(
                angle=45, color="slategrey", space_color="none"
            ),
            "color": "slategrey",
            "weight": 0.5,
        },
        tooltip="Exclusion buffer",
    ).add_to(m)

    # zones of interest
    folium.GeoJson(
        zones.to_crs(4326),
        name="Area of interest",
        style_function=lambda feature: {
            "fillColor": "none",
            "color": "darkslategrey",
            "weight": 0.5,
            "dashArray": "5, 5",
        },
        tooltip="Area of interest",
    ).add_to(m)

    # injection point
    folium.GeoJson(
        injection_point.to_crs(4326),
        name="Gas grid injection point",
        tooltip="Gas grid injection point<br>Dublin Port",
        marker=folium.Marker(icon=folium.Icon(icon="anchor", prefix="fa")),
        style_function=lambda feature: {"markerColor": "red"},
    ).add_to(m)

    m.add_child(cm1).add_child(cm2)
    m.add_child(BindColormap(fg1, cm1)).add_child(BindColormap(fg2, cm2))
    folium.LayerControl(collapsed=False).add_to(m)
    GroupedLayerControl(
        groups={"Caverns": [fg1, fg2]}, collapsed=False
    ).add_to(m)

    MeasureControl().add_to(m)
    MousePosition().add_to(m)
    Fullscreen().add_to(m)

    return m
