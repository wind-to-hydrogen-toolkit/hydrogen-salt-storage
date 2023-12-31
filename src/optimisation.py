"""optimisation.py
Functions to optimise storage locations for a wind farm with reference wind
turbines

NREL 15 MW reference turbine: https://doi.org/10.2172/1570430
"""

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import integrate

from src import functions as fns

# NREL 15 MW reference turbine specifications
REF_DIAMETER = 248  # m
REF_HUB_HEIGHT = 149  # m
REF_RATED_POWER = 15  # MW
REF_V_CUT_IN = 4  # m s-1
REF_V_RATED = 11  # m s-1
REF_V_CUT_OUT = 25  # m s-1

# power curve data
pc = pd.DataFrame(
    {
        "wind_speed": range(REF_V_RATED + 2),
        "power_curve": (
            [0] * 4 + [499, 1424, 2732, 4469, 6643, 9459, 12975] + [15000] * 2
        ),
    }
)
# for wind speeds between cut-in and rated
pc["gradient"] = pc["power_curve"] - pc["power_curve"].shift(1)
pc["intercept"] = pc["power_curve"] - pc["gradient"] * pc["wind_speed"]
pc = pc[
    (pc["wind_speed"] > REF_V_CUT_IN) & (pc["wind_speed"] < REF_V_RATED + 1)
]


def ref_power_curve(v: float) -> float:
    """
    Power curve for the reference wind turbine.

    NREL 15 MW reference turbine: https://doi.org/10.2172/1570430

    Parameters
    ----------
    v : Wind speed [m s-1]

    Returns
    -------
    - Wind turbine power output at a given wind speed from the power curve [MW]
    """

    if v < REF_V_CUT_IN:
        power_curve = 0
    elif REF_V_CUT_IN <= v < REF_V_RATED:
        pc_vals = pc[pc["wind_speed"] == np.trunc(v) + 1]
        power_curve = (pc_vals["gradient"] * v + pc_vals["intercept"]).iloc[0]
    elif REF_V_RATED <= v <= REF_V_CUT_OUT:
        power_curve = 15000
    elif v > REF_V_CUT_OUT:
        power_curve = 0

    return power_curve / 1000


def read_weibull_data(
    data_path_weibull: str,
    data_path_wind_farms: str,
    dat_extent: gpd.GeoSeries,
    dat_crs: int = fns.CRS,
) -> pd.DataFrame:
    """
    Extract Weibull parameters of wind speeds

    Parameters
    ----------
    data_path_weibull : Path to the Weibull parameter data Zip file
    data_path_weibull : Path to the wind farm data Zip file
    dat_extent : Extent of the Kish Basin data
    dat_crs : EPSG CRS

    Returns
    -------
    - Dataframe of k and c values for each wind farm
    """

    weibull_df = {}

    for w in ["c", "k"]:
        # read Weibull parameter data
        weibull_df[w] = fns.read_shapefile_from_zip(
            data_path=data_path_weibull, endswith=f"{w}_ITM.shp"
        )

        # read wind farm data
        wind_farms = fns.constraint_wind_farm(
            data_path=data_path_wind_farms,
            dat_extent=dat_extent,
        )

        # combine Codling wind farm polygons
        wind_farms["Name_"] = wind_farms["Name"].str.split(expand=True)[0]

        # convert CRS and keep areas intersecting with wind farms
        weibull_df[w] = (
            weibull_df[w]
            .to_crs(dat_crs)
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


def weibull_probability_distribution(k: float, c: float, v: float) -> float:
    """
    Equation (1) of Dinh et al. (2023a).
    https://doi.org/10.1016/j.ijhydene.2023.01.016

    Parameters
    ----------
    k : Shape (Weibull distribution parameter)
    c : Scale (Weibull distribution parameter) [m s-1]
    v : Wind speed [m s-1]

    Returns
    -------
    - Weibull probability distribution function [s m-1]
    """

    return k / c * np.power((v / c), (k - 1)) * np.exp(-np.power((v / c), k))


def annual_energy_production(
    n_turbines: int,
    k: float,
    c: float,
    w_loss: float = 0.1,
) -> float:
    """
    Equation (3) of Dinh et al. (2023a).
    https://doi.org/10.1016/j.ijhydene.2023.01.016

    In the integration, both the limit and absolute error tolerance have
    been increased.

    Parameters
    ----------
    n_turbines : Number of wind turbines in wind farm
    k : Shape (Weibull distribution parameter)
    c : Scale (Weibull distribution parameter) [m s-1]
    w_loss : Wake loss

    Returns
    -------
    - Annual energy production of wind farm [MWh]
    - Integral [MW]
    - Absolute error [MW]
    """

    integration = integrate.quad(
        lambda v: (
            ref_power_curve(v=v)
            * weibull_probability_distribution(k=k, c=c, v=v)
        ),
        a=REF_V_CUT_IN,
        b=REF_V_CUT_OUT,
        limit=100,  # an upper bound on the number of subintervals used
        epsabs=1.49e-6,  # absolute error tolerance
    )

    aep = 365 * 24 * n_turbines * (1 - w_loss) * integration[0]

    return aep, integration[0], integration[1]


def annual_hydrogen_production(
    aep: float,
    e_elec: float = 0.05,
    eta_conv: float = 0.93,
    e_pcl: float = 0.003,
) -> float:
    """
    Equation (4) of Dinh et al. (2023a).
    https://doi.org/10.1016/j.ijhydene.2023.01.016
    Based on Dinh et al. (2021).
    https://doi.org/10.1016/j.ijhydene.2020.04.232

    Constant values are based on Table 3 of Dinh et al. (2021) for PEM
    electrolysers in 2030.

    Parameters
    ----------
    aep : Annual energy production of wind farm [MWh]
    e_elec : Electricity required to supply the electrolyser to produce 1 kg
      of hydrogen [MWh kg-1]
    eta_conv : Conversion efficiency of the electrolyser
    e_pcl : Electricity consumed by other parts of the hydrogen plant
      [MWh kg-1]

    Returns
    -------
    - Annual hydrogen production [kg]
    """

    ahp = aep / ((e_elec / eta_conv) + e_pcl)

    return ahp


def capex_pipeline(
    e_cap,
    p_rate: float = 0.0055,
    rho: float = 8,
    v: float = 15,
) -> float:
    """
    Capital expenditure (CAPEX) for the pipeline.
    See Equation (18) of Dinh et al. (2023a):
    https://doi.org/10.1016/j.ijhydene.2023.01.016
    and section 3.1 of Dinh et al. (2023b).

    The estimation of offshore pipeline costs is based on the onshore pipeline
    calculations of Baufumé et al. (2013), multiplied by a factor of two to
    reflect the estimated cost scaling of onshore to expected offshore costs,
    as suggested by the International Renewable Energy Agency.

    In the reference case, the resulting capital expenditure for transmission
    pipelines and associated recompression stations was represented by a
    second order polynomial function.

    The electrolyser production rate was based on the Siemens - Silyzer 300.

    Because the electrolyser plant is assumed to be operating at full capacity
    at all times, the CAPEX was calculated considering a 75% utilization
    rate, i.e. the pipeline capacity is 33% oversized (International Energy
    Agency, 2019).

    - https://doi.org/10.1016/j.ijhydene.2012.12.147
    - https://www.irena.org/publications/2022/Apr/Global-hydrogen-trade-Part-II
    - https://www.iea.org/reports/the-future-of-hydrogen
    - https://assets.siemens-energy.com/siemens/assets/api/uuid:a193b68f-7ab4-4536-abe2-c23e01d0b526/datasheet-silyzer300.pdf

    Parameters
    ----------
    e_cap : Electrolyser capacity [MW]
    p_rate : Electrolyser production rate [kg s-1 MW-1]
    rho : Mass density of hydrogen [kg m-3]
    v : Average fluid velocity [m s-1]

    Returns
    -------
    - CAPEX of the pipeline per km of pipeline [€ km-1]
    """

    return (
        16000 * (e_cap * p_rate / (rho * v * np.pi)) +
        1197.2 * np.sqrt(e_cap * p_rate / (rho * v * np.pi)) +
        329
    ) * 1000 * 2


# def rotor_area() -> float:
#     """
#     Reference wind turbine rotor swept area, from Dinh et al. (2023a).
#     https://doi.org/10.1016/j.ijhydene.2023.01.016

#     Returns
#     -------
#     - Wind turbine rotor swept area [m2]
#     """

#     return np.pi * np.square(REF_DIAMETER) / 4


# def power_wind_resource(v: float, rho: float = 1.225) -> float:
#     """
#     Total power of the wind resource passing through the reference wind
#     turbine rotor

#     Parameters
#     ----------
#     v : Wind speed [m s-1]
#     rho : Air density [kg m-3]

#     Returns
#     -------
#     - Power contained in the wind resource [MW]
#     """

#     return 0.5 * rho * rotor_area() * np.power(v, 3) / 1000


# def power_coefficient(v: float) -> float:
#     """
#     Power coefficient curve of the reference wind turbine

#     Parameters
#     ----------
#     v : Wind speed [m s-1]

#     Returns
#     -------
#     - Power coefficient
#     """

#     try:
#         coeff = ref_power_curve(v=v) / power_wind_resource(v=v)
#     except ZeroDivisionError:
#         coeff = 0

#     return coeff


# def power_output_wind_turbine(v: float) -> float:
#     """
#     Equation (2) of Dinh et al. (2023a).
#     https://doi.org/10.1016/j.ijhydene.2023.01.016

#     Parameters
#     ----------
#     v : Wind speed [m s-1]

#     Returns
#     -------
#     - Power output of wind turbine [MW]
#     """

#     if v < REF_V_CUT_IN:
#         power_wt = 0
#     elif REF_V_CUT_IN <= v < REF_V_RATED:
#         power_wt = power_wind_resource(v=v) * power_coefficient(v=v)
#     elif REF_V_RATED <= v <= REF_V_CUT_OUT:
#         power_wt = REF_RATED_POWER
#     elif v > REF_V_CUT_OUT:
#         power_wt = 0

#     return power_wt


# def power_curve_weibull(k: float, c: float, v: float) -> float:
#     """
#     Power curve with Weibull equation applied

#     Parameters
#     ----------
#     k : Shape (Weibull distribution parameter)
#     c : Scale (Weibull distribution parameter) [m s-1]
#     v : Wind speed [m s-1]

#     Returns
#     -------
#     - Power curve multiplied by the Weibull probability distribution function
#       [MW s m-1 = M kg m s-2]
#     """

#     return (
#         ref_power_curve(v=v) *
#         weibull_probability_distribution(k=k, c=c, v=v)
#     )
