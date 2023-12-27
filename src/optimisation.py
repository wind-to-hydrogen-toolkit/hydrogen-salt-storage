"""optimisation.py
Functions to optimise storage locations
"""

import numpy as np
from scipy import integrate


def weibull_probability_distribution(k: float, c: float, v: float):
    """
    Equation (1) of Dinh et al. (2023).
    https://doi.org/10.1016/j.ijhydene.2023.01.016
    """

    return k / c * np.power((v / c), (k - 1)) * np.exp(-np.power((v / c), k))


def rotor_area(diameter: float):
    """
    Wind turbine rotor swept area, from Dinh et al. (2023).
    https://doi.org/10.1016/j.ijhydene.2023.01.016
    """

    return np.pi * np.square(diameter) / 4


def power_output_wind_turbine(
    v: float, area: float, c_power: float, v_i: float = 4, v_r: float = 11,
    v_o: float = 25, p_rated: float = 15
):
    """
    Equation (2) of Dinh et al. (2023).
    https://doi.org/10.1016/j.ijhydene.2023.01.016
    """

    if v < v_i:
        power_wt = 0
    elif v_i <= v < v_r:
        power_wt = 0.5 * 1.225 * area * c_power * np.power(v, 3)
    elif v_r <= v < v_o:
        power_wt = p_rated
    elif v >= v_o:
        power_wt = 0

    return power_wt


def annual_energy_production(
    n_turbines: int, diameter: float, c_power: float, k: float,
    c: float, v_i: float = 4, v_o: float = 25, w_loss: float = 0.1
):
    """
    Equation (3) of Dinh et al. (2023).
    https://doi.org/10.1016/j.ijhydene.2023.01.016
    """

    aep = (
        365 * 24 * n_turbines * (1 - w_loss) *
        integrate.quad(
            lambda v: (
                power_output_wind_turbine(
                    v=v, area=rotor_area(diameter=diameter), c_power=c_power
                ) *
                weibull_probability_distribution(k=k, c=c, v=v)
            ),
            v_i,
            v_o
        )
    )

    return aep
