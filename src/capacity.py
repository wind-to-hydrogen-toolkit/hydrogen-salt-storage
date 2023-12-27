"""capacity.py
Functions to calculate salt cavern volumes and storage capacities
"""

import numpy as np
import pandas as pd
from pyfluids import Fluid, FluidsList, Input


def cavern_volume(
    height: float, diameter: float = 80, theta: float = 20
) -> float:
    """
    Calculate the cavern volume. See Williams et al. (2022), eq. (1) and
    Jannel and Torquet (2022).

    - v_cavern = v_cylinder + 2 * v_cone
    - v_cylinder = pi * r^2 * h_cylinder = pi * r^2 * (h_cavern - 2 * h_cone)
    - v_cone = pi * r^2 * h_cone / 3 = pi * r^2 * (r * tan(theta)) / 3
    - v_cavern = pi * r^2 * (h_cavern - 4 / 3 * r * tan(theta))
    - v_cavern_corrected = v_cavern * scf * (1 - if * insf * bf)

    Parameters
    ----------
    height : Cavern height [m]
    diameter : Cavern diameter [m]
    theta : Cavern roof angle [deg]

    Returns
    -------
    - Corrected cavern volume [m3]
    """

    # calculate ideal cavern volume
    r = diameter / 2
    v_cavern = (
        np.pi * np.square(r) * (height - 4 / 3 * r * np.tan(np.deg2rad(theta)))
    )

    # apply correction factors
    v_cavern = v_cavern * 0.7 * (1 - 0.25 * 0.865 * 1.46)

    return v_cavern


def temperature_cavern_mid_point(
    height: float, depth_top: float, t_0: float = 10, t_delta: float = 37.5
) -> float:
    """
    Mid-point temperature. See Williams et al. (2022), eq. (2).

    t_delta value of 37.5 based on the geothermal gradient of Kish Basin wells
    in the Mercia Mudstone Group reported in English et al. (2023).

    Parameters
    ----------
    height : Cavern height [m]
    t_0 : Mean annual surface temperature [deg C]
    t_delta : Geothermal gradient; change in temperature with depth
        [deg C / km]
    depth_top : Cavern top depth [m]

    Returns
    -------
    - Mid-point temperature [K]
    """

    return t_0 + t_delta / 1000 * (depth_top + height / 2) + 273.15


def pressure_operating(
    thickness_overburden,
    rho_overburden: float = 2400,
    rho_salt: float = 2200,
    minf: float = 0.3,
    maxf: float = 0.8,
) -> tuple[float, float]:
    """
    Operating pressures. See Williams et al. (2022), eq. (3) and (4).

    Parameters
    ----------
    thickness_overburden : Overburden thickness / halite top depth [m]
    rho_overburden : Density of the overburden [kg m-3]
    rho_salt : Density of the halite [kg m-3]
    minf : Factor of lithostatic pressure for the minimum operating pressure
    maxf : Factor of lithostatic pressure for the maximum operating pressure

    Returns
    -------
    - Operating pressures [Pa]
    """

    # lithostatic pressure at the casing shoe
    # thickness of the overburden is the same as the depth to top of salt
    # thickness of salt above casing shoe = 50 m
    # acceleration due to gravity = 9.81
    p_casing = (rho_overburden * thickness_overburden + rho_salt * 50) * 9.81

    p_operating_min = minf * p_casing
    p_operating_max = maxf * p_casing

    return p_operating_min, p_operating_max


def density_hydrogen_gas(
    p_operating_min: float,
    p_operating_max: float,
    t_mid_point: float,
) -> tuple[float, float]:
    """
    Density of hydrogen at cavern conditions. See Williams et al. (2022),
    s. 3.4.2 and Caglayan et al. (2020), eq. (3).

    !Note: Check density values!

    http://www.coolprop.org/fluid_properties/fluids/Hydrogen.html

    rho = p * m / (z * r * t)
    where p: pressure [Pa], m: molar mass of hydrogen gas [kg mol-1], z:
    compressibility factor, r: universal gas constant [J K-1 mol-1], t:
    temperature [K]

    Parameters
    ----------
    p_operating_min : Minimum operating pressure [Pa]
    p_operating_max : Maximum operating pressure [Pa]
    t_mid_point : Mid-point temperature [K]

    Returns
    -------
    - Hydrogen gas density [kg m-3]
    """

    rho_h2 = []

    for p_min, p_max, t in zip(p_operating_min, p_operating_max, t_mid_point):
        h2_min = Fluid(FluidsList.Hydrogen).with_state(
            Input.pressure(p_min), Input.temperature(t)
        )

        h2_max = Fluid(FluidsList.Hydrogen).with_state(
            Input.pressure(p_max), Input.temperature(t)
        )

        rho_h2.append((h2_min.density, h2_max.density))

    # # rho / p = m / (z * r * t)
    # rho_p = (2.01588 / 1000) / (
    #     compressibility_factor * 8.314 * (t_mid_point)
    # )

    # rho_h2_min = p_operating_min * rho_p
    # rho_h2_max = p_operating_max * rho_p

    # rho_approx_min = (
    #     (2.01588 / 1000)
    #     / (h2_min.compressibility * 8.314 * (20 + 273.15))
    #     * 100e3
    # )

    rho_h2 = pd.DataFrame(rho_h2)

    return rho_h2[0], rho_h2[1]


def mass_hydrogen_working(
    rho_h2_min: float, rho_h2_max: float, v_cavern: float
) -> float:
    """
    The working mass of hydrogen in kg. See Williams et al. (2022), eq. (5)
    and (6).

    Parameters
    ----------
    rho_h2_min : Minimum hydrogen density [kg m-3]
    rho_h2_max : Maximum hydrogen density [kg m-3]
    v_cavern : Cavern volume [m3]

    Returns
    -------
    - Working mass of hydrogen [kg]
    """

    # stored mass of hydrogen at min and max operating pressures
    m_min_operating = rho_h2_min * v_cavern
    m_max_operating = rho_h2_max * v_cavern

    return m_max_operating - m_min_operating


def energy_storage_capacity(m_working: float, lhv: float = 119.96) -> float:
    """
    Cavern energy storage capacity in GWh. See Williams et al. (2022), eq. (7).

    Parameters
    ----------
    m_working : Working mass of hydrogen [kg]
    lhv : Lower heating value of hydrogen [MJ kg-1]

    Returns
    -------
    - Energy storage capacity of cavern [GWh]
    """

    return m_working * lhv / 3.6e6
