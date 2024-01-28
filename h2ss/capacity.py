"""Functions to calculate salt cavern volumes and storage capacities.

References
----------
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
        Ideal cavern volume [m³]

    Notes
    -----
    The cavern's ideal shape is a cylinder of height :math:`h_{cylinder}` [m]
    with two conical ends each with a height of :math:`h_{cone}` [m] based on
    [#Jannel22]_. Since the volume of a cylinder is defined as
    :math:`\\pi r^2 h_{cylinder}` and the volume of a cone is defined as
    :math:`\\frac{\\pi r^2 h_{cone}}{3}` (i.e. one third of a cylinder's
    volume), the bulk cavern volume, :math:`V_{bulk}` [m³] can therefore be
    deduced by calculating the volume of a cylinder of the cavern height,
    :math:`h` [m] and deducting 4/3 of the volume of a cone of height
    :math:`h_{cone}`. :math:`r` [m] is the radius of the cavern and
    :math:`h_{cone}` is therefore equivalent to :math:`r \\tan(\\theta)`,
    where :math:`\\theta` [°] is the cavern roof angle.

    .. math::
        V_{bulk} = \\pi \\cdot r^2 \\cdot h
        - \\frac{4}{3} \\, \\pi \\cdot r^2 \\cdot h_{cone}
    .. math::
        V_{bulk} = \\pi \\cdot r^2
        \\left(h - \\frac{4}{3} \\, h_{cone}\\right)
    .. math::
        V_{bulk} = \\pi \\cdot r^2 \\left(h - \\frac{4}{3} \\,
        r \\cdot \\tan(\\theta)\\right)
    """
    r = diameter / 2
    v_cavern = (
        np.pi * np.square(r) * (height - 4 / 3 * r * np.tan(np.deg2rad(theta)))
    )
    return v_cavern


def corrected_cavern_volume(
    v_cavern, f_scf=0.7, f_if=0.25, f_insf=0.865, f_bf=1.46
):
    """Apply correction factors to the cavern volume.

    Parameters
    ----------
    v_cavern : float
        Ideal cavern volume [m³]
    f_scf : float
        Shape correction factor
    f_if : float
        Insoluble fraction
    f_insf : float
        Correction for the fraction of insoluble material that remains in the
        cavern after mechanical sweeping
    f_bf : float
        Bulking factor

    Returns
    -------
    float
        Corrected cavern volume [m³]

    The corrected cavern volume, :math:`V_{cavern}` [m³] is approximated by
    applying several correction factors to the bulk cavern volume,
    :math:`V_{bulk}` [m³] as detailed in [#Williams22]_, Eqn. (1).

    .. math::
        V_{cavern} = V_{bulk} \\times SCF \\times (1 - IF \\times INSF
        \\times BF)
    .. math::
        V_{cavern} = V_{bulk} \\times 0.7 \\times (1 - 0.25 \\times 0.865
        \\times 1.46)
    .. math::
        V_{cavern} \\approx 0.48 \\, V_{bulk}

    where :math:`SCF` is the shape correction factor (deviation of the
    cavern's shape due to geological differences), :math:`IF` is the
    insoluble fraction (insoluble material within the salt), :math:`INSF`
    is the correction for the fraction of insoluble material that remains in
    the cavern after mechanical sweeping that forms the sump of the cavern,
    and :math:`BF` is the bulking factor (uneven stacking of the insoluble
    material that forms the cavern sump).
    """
    correction_factors = f_scf * (1 - f_if * f_insf * f_bf)
    return v_cavern * correction_factors


def temperature_cavern_mid_point(height, depth_top, t_0=10, delta_t=37.5):
    """Cavern mid-point temperature.

    Parameters
    ----------
    height : float
        Cavern height [m]
    depth_top : float
        Cavern top depth [m]
    t_0 : float
        Mean annual surface temperature [°C]
    delta_t : float
        Geothermal gradient; change in temperature with depth [°C km⁻¹]

    Returns
    -------
    float
        Mid-point temperature [K]

    Notes
    -----
    See [#Williams22]_, Eqn. (2).

    .. math::
        T_{midpoint} = T_0 + \\Delta_T \\,
        \\frac{(z + 0.5 \\, h)}{1,000} + 273.15

    where :math:`T_{midpoint}` is the cavern mid-point temperature [K],
    :math:`T_0` is the mean annual surface temperature [°C], :math:`\\Delta_T`
    is the change in temperature with depth, i.e. geothermal gradient
    [°C km⁻¹], :math:`z` is the cavern top depth [m], and :math:`h` is the
    cavern height [m].

    A :math:`\\Delta_T` of 37.5 °C km⁻¹ is used based on the geothermal
    gradient of Kish Basin wells in the Mercia Mudstone Group reported in
    [#English23]_.
    Note that the cavern top depth is used instead of the casing shoe depth,
    as the casing shoe is 30 m above the cavern top depth in this study,
    rather than at the cavern top as modelled in [#Williams22]_.
    """
    return t_0 + delta_t * (depth_top + height / 2) / 1000 + 273.15


def pressure_operating(
    thickness_overburden,
    thickness_salt=50,
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
    thickness_salt : float
        Thickness of the salt above the casing shoe [m]
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

    .. math::
        p_{casing} = (\\rho_{overburden} \\cdot t_{overburden} + \\rho_{salt}
        \\cdot t_{salt}) \\, g
    .. math::
        p_{H2min} = 0.3 \\, p_{casing}
    .. math::
        p_{H2max} = 0.8 \\, p_{casing}

    :math:`p_{casing}` is the lithostatic pressure at the casing shoe [Pa].
    The thickness of the overburden, :math:`t_{overburden}` [m] is the same
    as the depth to top of salt. :math:`t_{salt}` [m] is the thickness of
    the salt above the casing shoe. :math:`\\rho_{overburden}` and
    :math:`\\rho_{salt}` are the densities of the overburden and salt,
    respectively [kg m⁻³]. :math:`g` is the acceleration due to gravity
    [m s⁻²]. :math:`p_{H2min}` and :math:`p_{H2max}` are the minimum and
    maximum cavern operating pressures, respectively [Pa].
    """
    p_casing = (
        rho_overburden * thickness_overburden + rho_salt * thickness_salt
    ) * 9.81
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
        Min and max hydrogen gas densities [kg m⁻³]

    Notes
    -----
    This function uses the CoolProp [#Bell14]_ wrapper called PyFluids
    [#PyFluids]_. See [#Williams22]_, Section 3.4.2 and also the CoolProp
    documentation for useful information on hydrogen [#CoolProp]_.
    The ``pyproject.toml`` configuration file has been set such that the
    default units used by PyFluids are SI units. PyFluids can also be used to
    derive the compressibility factor.

    Based on [#Caglayan20]_, Eqn. (3), the density can be approximated using
    the following equation.

    .. math::
        \\rho = \\frac{p \\cdot M}{Z \\cdot R \\cdot T}

    where :math:`M` is the molar mass of hydrogen gas [kg mol⁻¹], :math:`R` is
    the universal gas constant [J K⁻¹ mol⁻¹], :math:`p` is the pressure [Pa],
    :math:`T` is the temperature [K], :math:`Z` is the compressibility factor,
    and :math:`\\rho` is the hydrogen gas density [kg m⁻³].
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
    rho_h2 = pd.DataFrame(rho_h2)
    return rho_h2[0], rho_h2[1]


def mass_hydrogen_working(rho_h2_min, rho_h2_max, v_cavern):
    """The working mass of hydrogen.

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
    tuple[float, float, float]
        Working mass and the minimum and maximum operating mass of hydrogen
        [kg]

    Notes
    -----
    See [#Williams22]_, Eqn. (5) and (6).
    The mass of hydrogen at the minimum operating pressure represents the
    cushion gas requirement.

    .. math::
        m_{min} = \\rho_{H2min} \\cdot V_{cavern}
    .. math::
        m_{max} = \\rho_{H2max} \\cdot V_{cavern}
    .. math::
        m_{working} = m_{max} - m_{min}

    The working mass of hydrogen, :math:`m_{working}` [kg] is the difference
    between the stored mass of hydrogen at maximum, :math:`m_{max}` and
    minimum, :math:`m_{min}` operating pressures [kg], which were derived
    using the minimum, :math:`\\rho_{min}` and maximum, :math:`\\rho_{max}`
    hydrogen densities [kg m⁻³], respectively, and the cavern volume
    :math:`V_{cavern}`, [m³].
    """
    m_min_operating = rho_h2_min * v_cavern
    m_max_operating = rho_h2_max * v_cavern
    m_working = m_max_operating - m_min_operating
    return m_working, m_min_operating, m_max_operating


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
    The energy storage capacity of the cavern, :math:`E_{cavern}` [GWh] is a
    function of the working mass of hydrogen, :math:`m_{working}` [kg] and the
    lower heating value of hydrogen, :math:`LHV` [MJ kg⁻¹].

    .. math::
        E_{cavern} = m_{working} \\, \\frac{LHV}{3,600,000}
    """
    return m_working * lhv / 3.6e6


def calculate_capacity_dataframe(cavern_df):
    """Calculate volumes and storage capacities for a dataframe of caverns.

    Parameters
    ----------
    cavern_df : gpd.GeoDataFrame
        Dataframe of caverns

    Returns
    -------
    gpd.GeoDataFrame
        Updated dataframe of caverns with volume and capacity values
    """
    cavern_df["cavern_volume"] = cavern_volume(
        height=cavern_df["cavern_height"]
    )
    cavern_df["cavern_volume"] = corrected_cavern_volume(
        v_cavern=cavern_df["cavern_volume"]
    )

    cavern_df["t_mid_point"] = temperature_cavern_mid_point(
        height=cavern_df["cavern_height"], depth_top=cavern_df["cavern_depth"]
    )

    (
        cavern_df["p_operating_min"],
        cavern_df["p_operating_max"],
    ) = pressure_operating(thickness_overburden=cavern_df["TopDepthSeabed"])

    cavern_df["rho_min"], cavern_df["rho_max"] = density_hydrogen_gas(
        p_operating_min=cavern_df["p_operating_min"],
        p_operating_max=cavern_df["p_operating_max"],
        t_mid_point=cavern_df["t_mid_point"],
    )

    (
        cavern_df["working_mass"],
        cavern_df["mass_operating_min"],
        cavern_df["mass_operating_max"],
    ) = mass_hydrogen_working(
        rho_h2_min=cavern_df["rho_min"],
        rho_h2_max=cavern_df["rho_max"],
        v_cavern=cavern_df["cavern_volume"],
    )

    cavern_df["capacity"] = energy_storage_capacity(
        m_working=cavern_df["working_mass"]
    )

    return cavern_df
