"""optimisation.py
Functions to optimise storage locations for a wind farm with reference wind
turbines

NREL 15 MW reference turbine: https://doi.org/10.2172/1570430
"""

import numpy as np
import pandas as pd
from scipy import integrate

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
        "wind_speed": range(REF_V_CUT_OUT + 1),
        "power_curve": (
            [0] * 4 + [499, 1424, 2732, 4469, 6643, 9459, 12975] + [15000] * 15
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


# def rotor_area() -> float:
#     """
#     Reference wind turbine rotor swept area, from Dinh et al. (2023).
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


def weibull_probability_distribution(k: float, c: float, v: float) -> float:
    """
    Equation (1) of Dinh et al. (2023).
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


# def power_output_wind_turbine(v: float) -> float:
#     """
#     Equation (2) of Dinh et al. (2023).
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
#       [MW s m-1 / M kg m s-2]
#     """

#     return (
#         ref_power_curve(v=v) *
#         weibull_probability_distribution(k=k, c=c, v=v)
#     )


def annual_energy_production(
    n_turbines: int,
    k: float,
    c: float,
    w_loss: float = 0.1,
) -> float:
    """
    Equation (3) of Dinh et al. (2023).
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
    - Annual energy production of wind farm [GWh]
    """

    integration = integrate.quad(
        lambda v: (
            ref_power_curve(v=v)
            * weibull_probability_distribution(k=k, c=c, v=v)
        ),
        a=REF_V_CUT_IN,
        b=REF_V_CUT_OUT,
        limit=100,  # an upper bound on the number of subintervals used
        epsabs=1.49e-6  # absolute error tolerance
    )

    aep = 365 * 24 * n_turbines * (1 - w_loss) * integration[0]

    return aep / 1000, integration


def annual_hydrogen_production(
    n_turbines: int,
    k: float,
    c: float,
    e_elec: float,
    eta_conv: float,
    e_pcl: float,
) -> float:
    """
    Equation (4) of Dinh et al. (2023).
    https://doi.org/10.1016/j.ijhydene.2023.01.016

    Parameters
    ----------
    n_turbines : Number of wind turbines in wind farm
    k : Shape (Weibull distribution parameter)
    c : Scale (Weibull distribution parameter) [m s-1]
    e_elec : Electricity required to supply the electrolyser to produce 1 kg
      of hydrogen [MWh/kg]
    eta_conv : Conversion efficiency of the electrolyser
    e_pcl : Electricity consumed by other parts of the hydrogen plant [MWh/kg]

    Returns
    -------
    - Annual hydrogen production [kg]
    """

    ahp = annual_energy_production(n_turbines=n_turbines, k=k, c=c) / (
        (e_elec / eta_conv) + e_pcl
    )

    return ahp
