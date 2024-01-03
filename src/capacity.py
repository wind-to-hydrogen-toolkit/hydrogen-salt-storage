"""Functions to calculate salt cavern volumes and storage capacities

.. rubric:: References
.. [#Jannel22] Jannel, H. and Torquet, M. (2022). Conceptual design of salt
    cavern and porous media underground storage site. Hystories deliverable
    D7.1-1. Hystories. Available at:
    https://hystories.eu/wp-content/uploads/2022/05/Hystories_D7.1-1-Conceptual-design-of-salt-cavern-and-porous-media-underground-storage-site.pdf
    (Accessed: 9 October 2023).
.. [#Williams22] Williams, J. D. O., Williamson, J. P., Parkes, D., Evans, D.
    J., Kirk, K. L., Sunny, N., Hough, E., Vosper, H., and Akhurst, M. C.
    (2022). ‘Does the United Kingdom have sufficient geological storage
    capacity to support a hydrogen economy? Estimating the salt cavern storage
    potential of bedded halite formations’, Journal of Energy Storage, 53, p.
    105109. https://doi.org/10.1016/j.est.2022.105109.
.. [#English23] English, J. M., English, K. L., Dunphy, R. B., Blake, S.,
    Walsh, J., Raine, R., Vafeas, N. A., and Salgado, P. R. (2023). ‘An
    Overview of Deep Geothermal Energy and Its Potential on the Island of
    Ireland’, First Break, 41(2), pp. 33–43.
    https://doi.org/10.3997/1365-2397.fb2023009.
.. [#Caglayan20] Caglayan, D. G., Weber, N., Heinrichs, H. U., Linßen, J.,
    Robinius, M., Kukla, P. A., and Stolten, D. (2020). ‘Technical potential
    of salt caverns for hydrogen storage in Europe’, International Journal of
    Hydrogen Energy, 45(11), pp. 6793–6805.
    https://doi.org/10.1016/j.ijhydene.2019.12.161.
.. [#Bell14] Bell, I. H., Wronski, J., Quoilin, S., and Lemort, V. (2014).
    ‘Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the
    Open-Source Thermophysical Property Library CoolProp’, Industrial &
    Engineering Chemistry Research, 53(6), pp. 2498–2508.
    https://doi.org/10.1021/ie4033999.
.. [#PyFluids] Portyanikhin, V. (2023). ‘portyanikhin/PyFluids’. Available at:
    https://github.com/portyanikhin/PyFluids (Accessed: 1 January 2024).
.. [#CoolProp] Bell, I. H. and the CoolProp Team (n.d.). Hydrogen - CoolProp
    documentation. Available at:
    http://www.coolprop.org/fluid_properties/fluids/Hydrogen.html
    (Accessed: 10 December 2023).
"""

import numpy as np
import pandas as pd
from pyfluids import Fluid, FluidsList, Input


def cavern_volume(height, diameter=80, theta=20):
    """Calculate the cavern volume.

    Parameters
    ----------
    height : float
        Cavern height [m]
    diameter : float
        Cavern diameter [m]
    theta : float
        Cavern roof angle [°]

    Returns
    -------
    float
        Corrected cavern volume [m³]

    Notes
    -----
    The cavern's ideal shape is a cylinder of height :math:`h_{cylinder}` [m]
    with two conical ends each with a height of :math:`h_{cone}` [m] based on
    [#Jannel22]_. Since the volume of a cylinder is defined as
    :math:`\\pi r^2 h` and the volume of a cone is defined as
    :math:`\\frac{\\pi r^2 h}{3}` (i.e. one third of a cylinder's volume), the
    ideal cavern volume, :math:`V_{ideal}` [m³] can therefore be derived by
    calculating the volume of a cylinder of the cavern height,
    :math:`h_{cavern}` [m] and deducting 4/3 of the volume of a cone of height
    :math:`h_{cone}`. :math:`r` [m] is the radius of the cavern and
    :math:`h_{cone}` is therefore equivalent to :math:`r \\tan(\\theta)`,
    where :math:`\\theta` [°] is the cavern roof angle.

    .. math::
        V_{ideal} = \\pi \\cdot r^2 \\cdot h_{cavern}
        - \\frac{4}{3} \\times \\pi \\cdot r^2 \\cdot h_{cone}
    .. math::
        V_{ideal} = \\pi \\cdot r^2
        \\left(h_{cavern} - \\frac{4}{3} \\times h_{cone}\\right)
    .. math::
        V_{ideal} = \\pi \\cdot r^2 \\left(h_{cavern} - \\frac{4}{3} \\times
        r \\cdot \\tan(\\theta)\\right)

    The corrected cavern volume, :math:`V_{cavern}` [m³] is approximated by
    applying several correction factors as detailed in [#Williams22]_, Eqn.
    (1).

    .. math::
        V_{cavern} = V_{ideal} \\times SCF \\times (1 - IF \\times INSF
        \\times BF)
    .. math::
        V_{cavern} = V_{ideal} \\times 0.7 \\times (1 - 0.25 \\times 0.865
        \\times 1.46)
    .. math::
        V_{cavern} \\approx 0.48 \\times V_{ideal}
    """
    # calculate ideal cavern volume
    r = diameter / 2
    v_cavern = (
        np.pi * np.square(r) * (height - 4 / 3 * r * np.tan(np.deg2rad(theta)))
    )

    # apply correction factors
    correction_factors = 0.7 * (1 - 0.25 * 0.865 * 1.46)
    v_cavern = v_cavern * correction_factors

    return v_cavern


def temperature_cavern_mid_point(height, depth_top, t_0=10, t_delta=37.5):
    """Cavern mid-point temperature.

    Parameters
    ----------
    height : float
        Cavern height [m]
    t_0 : float
        Mean annual surface temperature [°C]
    t_delta : float
        Geothermal gradient; change in temperature with depth [°C km⁻¹]
    depth_top : float
        Cavern top depth [m]

    Returns
    -------
    float
        Mid-point temperature [K]

    Notes
    -----
    See [#Williams22]_, Eqn. (2).

    A :math:`T_\\Delta` [°C km⁻¹] value of 37.5 is used based on the geothermal
    gradient of Kish Basin wells in the Mercia Mudstone Group reported in
    [#English23]_.
    Note that the cavern top depth is used instead of the casing shoe depth,
    as the casing shoe is 30 m above the cavern top depth in this study,
    rather than at the cavern top.
    """
    return t_0 + t_delta * (depth_top + height / 2) / 1000 + 273.15


def pressure_operating(
    thickness_overburden,
    rho_overburden=2400,
    rho_salt=2200,
    minf=0.3,
    maxf=0.8,
):
    """Cavern operating pressures.

    Parameters
    ----------
    thickness_overburden : float
        Overburden thickness / halite top depth [m]
    rho_overburden : float
        Density of the overburden [kg m⁻³]
    rho_salt : float
        Density of the halite [kg m⁻³]
    minf : float
        Factor of lithostatic pressure for the minimum operating pressure
    maxf : float
        Factor of lithostatic pressure for the maximum operating pressure

    Returns
    -------
    tuple[float, float]
        Min and max operating pressures [Pa]

    Notes
    -----
    See [#Williams22]_, Eqn. (3) and (4).
    """
    # lithostatic pressure at the casing shoe
    # thickness of the overburden is the same as the depth to top of salt
    # thickness of salt above casing shoe = 50 m
    # acceleration due to gravity = 9.81
    p_casing = (rho_overburden * thickness_overburden + rho_salt * 50) * 9.81

    p_operating_min = minf * p_casing
    p_operating_max = maxf * p_casing

    return p_operating_min, p_operating_max


def density_hydrogen_gas(p_operating_min, p_operating_max, t_mid_point):
    """Density of hydrogen at cavern conditions.

    Parameters
    ----------
    p_operating_min : float
        Minimum operating pressure [Pa]
    p_operating_max : float
        Maximum operating pressure [Pa]
    t_mid_point : float
        Mid-point temperature [K]

    Returns
    -------
    tuple[float, float]
        Min and max hydrogen gas density [kg m⁻³]

    Notes
    -----
    This function uses the CoolProp [#Bell14]_ wrapper called PyFluids
    [#PyFluids]_. See [#Williams22]_, Section 3.4.2 and also the CoolProp
    documentation for useful information on hydrogen [#CoolProp]_.
    The `pyfluids.ini` configuration file has been set such that the default
    units used by PyFluids are SI units. PyFluids can also be used to derive
    the compressibility factor.

    Based on [#Caglayan20]_, Eqn. (3), the density can be approximated using
    the following equation.

    .. math::
        \\rho = \\frac{P \\cdot M}{Z \\cdot R \\cdot T}

    where :math:`M` is the molar mass of hydrogen gas [kg mol⁻¹], :math:`R` is
    the universal gas constant [J K⁻¹ mol⁻¹], :math:`P` is the pressure [Pa],
    :math:`T` is the temperature [K], :math:`Z` is the compressibility factor,
    and :math:`\\rho` is the hydrogen gas density [kg m⁻³].
    """
    rho_h2 = []

    for p_min, p_max, t in zip(p_operating_min, p_operating_max, t_mid_point):
        if p_min < p_max:
            h2_min = Fluid(FluidsList.Hydrogen).with_state(
                Input.pressure(p_min), Input.temperature(t)
            )

            h2_max = Fluid(FluidsList.Hydrogen).with_state(
                Input.pressure(p_max), Input.temperature(t)
            )

            rho_h2.append((h2_min.density, h2_max.density))
        else:
            print(
                "Error: p_operating_min seems to be larger than"
                "p_operating_max!"
            )

    rho_h2 = pd.DataFrame(rho_h2)

    return rho_h2[0], rho_h2[1]


def mass_hydrogen_working(rho_h2_min, rho_h2_max, v_cavern):
    """The working mass of hydrogen in kg.

    Parameters
    ----------
    rho_h2_min : float
        Minimum hydrogen density [kg m⁻³]
    rho_h2_max : float
        Maximum hydrogen density [kg m⁻³]
    v_cavern : float
        Cavern volume [m³]

    Returns
    -------
    float
        Working mass of hydrogen [kg]

    Notes
    -----
    See [#Williams22]_, Eqn. (5) and (6).
    """
    # stored mass of hydrogen at min and max operating pressures
    if rho_h2_min < rho_h2_max:
        m_min_operating = rho_h2_min * v_cavern
        m_max_operating = rho_h2_max * v_cavern
    else:
        print("Error: rho_h2_min seems to be larger than rho_h2_max!")

    return m_max_operating - m_min_operating


def energy_storage_capacity(m_working, lhv=119.96):
    """Cavern energy storage capacity.

    Parameters
    ----------
    m_working : float
        Working mass of hydrogen [kg]
    lhv : float
        Lower heating value of hydrogen [MJ kg⁻¹]

    Returns
    -------
    float
        Energy storage capacity of cavern [GWh]

    Notes
    -----
    See [#Williams22]_, Eqn. (7).
    The energy storage capacity of the cavern, :math:`E` [GWh] is a function
    of the working mass of hydrogen, :math:`m_{working}` [kg] and the lower
    heating value of hydrogen, :math:`LHV` [MJ kg⁻¹].

    .. math::
        E = m_{working} \\times \\frac{LHV}{3,600,000}
    """
    return m_working * lhv / 3.6e6
