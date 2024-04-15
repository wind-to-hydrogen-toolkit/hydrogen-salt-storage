"""Functions to compare and validate results.

References
----------
.. [#Deane21] Deane, P. (2021) Our Climate Neutral Future: Zero by 50. Wind
    Energy Ireland. Available at:
    https://windenergyireland.com/images/files/our-climate-neutral-future-0by50-final-report.pdf
    (Accessed: 8 February 2024).
.. [#Jannel22] Jannel, H. and Torquet, M. (2022). Conceptual design of salt
    cavern and porous media underground storage site. Hystories deliverable
    D7.1-1. Hystories. Available at:
    https://hystories.eu/wp-content/uploads/2022/05/Hystories_D7.1-1-Conceptual-design-of-salt-cavern-and-porous-media-underground-storage-site.pdf
    (Accessed: 9 October 2023).
"""

import os
import sys
from h2ss import capacity as cap
from h2ss import data as rd
from h2ss import functions as fns
# from h2ss import optimisation as opt


class HiddenPrints:
    """Suppress print statements: https://stackoverflow.com/a/45669280"""

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def electricity_demand_ie(caverns_df):
    """Compare the total capacity to Ireland's electricity demand in 2050.

    Parameters
    ----------
    cavern_df : geopandas.GeoDataFrame
        Geodataframe of caverns within the zone of interest

    Notes
    -----
    Figures from [#Deane21]_.
    """
    print(
        "Energy capacity as a percentage of Ireland's electricity demand "
        "in 2050:",
        f"{(caverns_df['capacity'].sum() / 1000 / 122 * 100):.2f}–"
        f"{(caverns_df['capacity'].sum() / 1000 / 84 * 100):.2f}%"
    )


def cavern_volumes(
    cavern_df, volume_case=380000, minimum_fraction=0.85
):
    """Verify whether cavern volumes are within recommended ranges.

    Parameters
    ----------
    cavern_df : geopandas.GeoDataFrame
        Geodataframe of caverns within the zone of interest
    volume_case : float
        Cavern volume corresponding to a Hystories Project investment scenario
    minimum_fraction : float
        The fraction of ``volume_case`` that is allowed as the minimum

    Returns
    -------
    geopandas.GeoDataFrame
        Dataframe of available caverns

    Notes
    -----
    See [#Jannel22]_ for the Hystories Project investment scenarios. The
    volume used is the free gas volume of a cavern. The volume should be no
    less than 85% of the case's volume.
    """
    return cavern_df[
        cavern_df["cavern_volume"] >= volume_case * minimum_fraction
    ]


def load_all_data():
    ds, extent = rd.kish_basin_data_depth_adjusted(
        dat_path=os.path.join("data", "kish-basin"),
        bathymetry_path=os.path.join("data", "bathymetry"),
    )

    exclusions = {}

    # exploration wells
    _, exclusions["wells_b"] = fns.constraint_exploration_well(
        data_path=os.path.join(
            "data",
            "exploration-wells",
            "Exploration_Wells_Irish_Offshore.shapezip.zip",
        )
    )

    # wind farms
    exclusions["wind_farms"] = fns.constraint_wind_farm(
        data_path=os.path.join(
            "data", "wind-farms", "wind-farms-foreshore-process.zip"
        ),
        dat_extent=extent,
    )

    # frequent shipping routes
    _, exclusions["shipping_b"] = fns.constraint_shipping_routes(
        data_path=os.path.join(
            "data", "shipping", "shipping_frequently_used_routes.zip"
        ),
        dat_extent=extent,
    )

    # shipwrecks
    _, exclusions["shipwrecks_b"] = fns.constraint_shipwrecks(
        data_path=os.path.join(
            "data", "shipwrecks", "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"
        ),
        dat_extent=extent,
    )

    # subsea cables
    _, exclusions["cables_b"] = fns.constraint_subsea_cables(
        data_path=os.path.join("data", "subsea-cables", "KIS-ORCA.gpkg")
    )

    return ds, extent, exclusions


def capacity_function(ds, extent, exclusions, cavern_diameter, min_cavern_height):
    """
    """

    # distance from salt formation edge
    edge_buffer = fns.constraint_halite_edge(dat_xr=ds, buffer=cavern_diameter * 3)

    zones, zds = fns.zones_of_interest(
        dat_xr=ds, constraints={"net_height": min_cavern_height, "min_depth": 500, "max_depth": 2000},
    )

    caverns = fns.generate_caverns_hexagonal_grid(
        zones_df=zones,
        dat_extent=extent,
        diameter=cavern_diameter,
        separation=cavern_diameter * 4,
    )

    caverns = fns.cavern_dataframe(
        dat_zone=zds,
        cavern_df=caverns,
        depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
    )

    # label caverns by depth and heights
    caverns = fns.label_caverns(
        cavern_df=caverns,
        heights=[min_cavern_height],
        depths={"min": 500, "min_opt": 1000, "max_opt": 1500, "max": 2000},
    )

    with HiddenPrints():
        caverns, _ = fns.generate_caverns_with_constraints(
            cavern_df=caverns,
            exclusions={
                "wells": exclusions["wells_b"],
                "wind_farms": exclusions["wind_farms"],
                "shipwrecks": exclusions["shipwrecks_b"],
                "shipping": exclusions["shipping_b"],
                "cables": exclusions["cables_b"],
                "edge": edge_buffer,
            },
        )

    caverns["cavern_total_volume"] = cap.cavern_volume(
        height=caverns["cavern_height"], diameter=cavern_diameter
    )
    caverns["cavern_volume"] = cap.corrected_cavern_volume(
        v_cavern=caverns["cavern_total_volume"], f_if=0
    )

    caverns["t_mid_point"] = cap.temperature_cavern_mid_point(
        height=caverns["cavern_height"], depth_top=caverns["cavern_depth"]
    )

    (
        caverns["p_operating_min"],
        caverns["p_operating_max"],
    ) = cap.pressure_operating(thickness_overburden=caverns["TopDepthSeabed"], depth_water=-caverns["Bathymetry"])

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

    caverns["cavern_diameter"] = cavern_diameter

    df = caverns[["cavern_diameter", "cavern_height", "capacity"]].copy()

    return df
