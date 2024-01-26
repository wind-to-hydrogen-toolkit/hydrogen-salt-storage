"""Functions to optimise hydrogen storage locations.

References
----------
.. [#Musial19] Musial, W. D., Beiter, P. C., Nunemaker, J., Heimiller, D. M.,
    Ahmann, J., and Busch, J. (2019). Oregon Offshore Wind Site Feasibility
    and Cost Study. Technical Report NREL/TP-5000-74597. Golden, CO: National
    Renewable Energy Laboratory. https://doi.org/10.2172/1570430.
.. [#Dinh23a] Dinh, Q. V., Dinh, V. N., Mosadeghi, H., Todesco Pereira, P. H.,
    and Leahy, P. G. (2023). ‘A geospatial method for estimating the
    levelised cost of hydrogen production from offshore wind’,
    International Journal of Hydrogen Energy, 48(40), pp. 15000–15013.
    https://doi.org/10.1016/j.ijhydene.2023.01.016.
.. [#Pryor18] Pryor, S. C., Shepherd, T. J., and Barthelmie, R. J. (2018).
    ‘Interannual variability of wind climates and wind turbine annual energy
    production’, Wind Energy Science, 3(2), pp. 651–665.
    https://doi.org/10.5194/wes-3-651-2018.
.. [#Dinh21] Dinh, V. N., Leahy, P., McKeogh, E., Murphy, J., and Cummins, V.
    (2021). ‘Development of a viability assessment model for hydrogen
    production from dedicated offshore wind farms’, International Journal of
    Hydrogen Energy, 46(48), pp. 24620–24631.
    https://doi.org/10.1016/j.ijhydene.2020.04.232.
.. [#Dinh23b] Dinh, Q. V., Dinh, V. N., and Leahy, P. (2023). ‘A method to map
    the levelised cost of hydrogen from offshore wind farms coupled to onshore
    electrolysers via HVDC’, IOP Conference Series: Earth and Environmental
    Science, 1281(1), p. 012005.
    https://doi.org/10.1088/1755-1315/1281/1/012005.
.. [#Dinh23c] Dinh, Q. V., Todesco Pereira, P. H., Dinh, V. N., Nagle, A. J.,
    and Leahy, P. G. (2023). ‘Optimising the levelised cost of transmission
    for green hydrogen and ammonia in new-build offshore energy
    infrastructure: pipelines, tankers, and HVDC’.
.. [#Baufume13] Baufumé, S., Grüger, F., Grube, T., Krieg, D., Linssen, J.,
    Weber, M., Hake, J.-F., and Stolten, D. (2013). ‘GIS-based scenario
    calculations for a nationwide German hydrogen pipeline infrastructure’,
    International Journal of Hydrogen Energy, 38(10), pp. 3813–3829.
    https://doi.org/10.1016/j.ijhydene.2012.12.147.
.. [#IRENA22] International Renewable Energy Agency (2022). Global Hydrogen
    Trade to Meet the 1.5°C Climate Goal: Technology Review of Hydrogen
    Carriers. Abu Dhabi: International Renewable Energy Agency.
    ISBN 978-92-9260-431-8. Available at:
    https://www.irena.org/publications/2022/Apr/Global-hydrogen-trade-Part-II
    (Accessed: 31 December 2023).
.. [#Siemens] Siemens Energy (n.d.). Silyzer 300. Datasheet. Erlangen.
    Available at:
    https://assets.siemens-energy.com/siemens/assets/api/uuid:a193b68f-7ab4-4536-abe2-c23e01d0b526/datasheet-silyzer300.pdf
    (Accessed: 31 December 2023).
.. [#IEA19] International Energy Agency (2019). The Future of Hydrogen:
    Seizing today’s opportunities. Paris: International Energy Agency.
    Available at: https://www.iea.org/reports/the-future-of-hydrogen
    (Accessed: 14 August 2023).
"""

from functools import partial

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
        "wind_speed": range(REF_V_RATED + 2),
        "power_curve": (
            [0] * 4
            + [0.499, 1.424, 2.732, 4.469, 6.643, 9.459, 12.975]
            + [REF_RATED_POWER] * 2
        ),
    }
)
# for wind speeds between cut-in and rated
pc["gradient"] = pc["power_curve"] - pc["power_curve"].shift(1)
pc["intercept"] = pc["power_curve"] - pc["gradient"] * pc["wind_speed"]
pc = pc[
    (pc["wind_speed"] > REF_V_CUT_IN) & (pc["wind_speed"] < REF_V_RATED + 1)
]


def ref_power_curve(v):
    """Power curve for the reference wind turbine.

    Parameters
    ----------
    v : float
        Wind speed [m s⁻¹]

    Returns
    -------
    float
        Wind turbine power output at a given wind speed from the power curve
        [MW]

    Notes
    -----
    The NREL 15 MW reference wind turbine is used; see [#Musial19]_, Appendix
    D, p. 84, for the data used to generate the power curve.
    """
    if v < REF_V_CUT_IN:
        power_curve = 0
    elif REF_V_CUT_IN <= v < REF_V_RATED:
        pc_vals = pc[pc["wind_speed"] == np.trunc(v) + 1]
        power_curve = (pc_vals["gradient"] * v + pc_vals["intercept"]).iloc[0]
    elif REF_V_RATED <= v <= REF_V_CUT_OUT:
        power_curve = REF_RATED_POWER
    elif v > REF_V_CUT_OUT:
        power_curve = 0
    return power_curve


def weibull_probability_distribution(v, k, c):
    """Weibull probability distribution function.

    Parameters
    ----------
    v : float
        Wind speed [m s⁻¹]
    k : float
        Shape (Weibull distribution parameter, k)
    c : float
        Scale (Weibull distribution parameter, C) [m s⁻¹]

    Returns
    -------
    float
        Weibull probability distribution function [s m⁻¹]

    Notes
    -----
    The Weibull probability distribution function, :math:`f(v)` [s m⁻¹] is
    based on Eqn. (1) of [#Dinh23a]_, where :math:`k` and :math:`C`
    [m s⁻¹] are the shape and scale Weibull distribution parameters,
    respectively, and :math:`v` is the wind speed.

    .. math::
        f(v) = \\frac{k}{C} \\, \\left( \\frac{v}{C} \\right)^{k-1}
        \\cdot \\exp \\left( -\\left( \\frac{v}{C} \\right)^k \\right)

    See also [#Pryor18]_, Eqn. (2).
    """
    return k / c * np.power((v / c), (k - 1)) * np.exp(-np.power((v / c), k))


def weibull_power_curve(v, k, c):
    """Weibull probability distribution function multiplied by the power curve.

    Parameters
    ----------
    v : float
        Wind speed between cut-in and cut-out [m s⁻¹]
    k : float
        Shape (Weibull distribution parameter, k)
    c : float
        Scale (Weibull distribution parameter, C) [m s⁻¹]

    Returns
    -------
    float
        Weibull probability distribution multiplied by the power curve
        [MW s m⁻¹]
    """
    if REF_V_CUT_IN <= v < REF_V_RATED:
        pc_vals = pc[pc["wind_speed"] == np.trunc(v) + 1]
        power_curve = (pc_vals["gradient"] * v + pc_vals["intercept"]).iloc[
            0
        ] * weibull_probability_distribution(v=v, k=k, c=c)
    elif REF_V_RATED <= v <= REF_V_CUT_OUT:
        power_curve = REF_RATED_POWER * weibull_probability_distribution(
            v=v, k=k, c=c
        )
    return power_curve


def number_of_turbines(owf_cap, wt_power=REF_RATED_POWER):
    """Number of reference wind turbines in the offshore wind farm.

    Parameters
    ----------
    owf_cap : float
        Maximum nameplate capacity of the proposed offshore wind farm [MW]
    wt_power : float
        Rated power of the reference wind turbine [MW]

    Returns
    -------
    int
        Number of wind turbines in the offshore wind farm comprising of
        reference wind turbines
    """
    return (owf_cap / wt_power).astype(int)


def annual_energy_production(n_turbines, k, c, w_loss=0.1):
    """Annual energy production of the wind farm.

    Parameters
    ----------
    n_turbines : int
        Number of wind turbines in wind farm
    k : float
        Shape (Weibull distribution parameter)
    c : float
        Scale (Weibull distribution parameter) [m s⁻¹]
    w_loss : float
        Wake loss

    Returns
    -------
    tuple[float, float, float]
        Annual energy production of wind farm [MWh], integral [MW], and
        absolute error [MW]

    Notes
    -----
    The annual energy production, :math:`AEP` [MWh], is based on Eqn. (3)
    of [#Dinh23a]_, where :math:`n` is the number of turbines in the wind
    farm, :math:`w` is the wake loss, which is assumed to be a constant
    value of 0.1, :math:`v_i` and :math:`v_o` [m s⁻¹] are the cut-in and
    cut-out speeds of the wind turbine, respectively, :math:`P(v)` [MW] is
    the wind turbine power output, and :math:`f(v)` [s m⁻¹] is the Weibull
    probability distribution function.

    .. math::
        AEP = 365 \\times 24 \\times n \\times
        \\left( 1 - w \\right) \\times
        \\int\\limits_{v_i}^{v_o} P(v) \\, f(v) \\,\\mathrm{d}v

    In the function's implementation, both the limit and absolute error
    tolerance for the integration have been increased.
    """
    integration = integrate.quad(
        partial(weibull_power_curve, k=k, c=c),
        a=REF_V_CUT_IN,
        b=REF_V_CUT_OUT,
        limit=100,  # an upper bound on the number of subintervals used
        epsabs=1.49e-6,  # absolute error tolerance
    )
    aep = 365 * 24 * n_turbines * (1 - w_loss) * integration[0]
    return aep, integration[0], integration[1]


def annual_hydrogen_production(aep, e_elec=0.05, eta_conv=0.93, e_pcl=0.003):
    """Annual hydrogen production from the wind farm's energy.

    Parameters
    ----------
    aep : float
        Annual energy production of wind farm [MWh]
    e_elec : float
        Electricity required to supply the electrolyser to produce 1 kg of
        hydrogen [MWh kg⁻¹]
    eta_conv : float
        Conversion efficiency of the electrolyser
    e_pcl : float
        Electricity consumed by other parts of the hydrogen plant [MWh kg⁻¹]

    Returns
    -------
    float
        Annual hydrogen production [kg]

    Notes
    -----
    Eqn. (4) of [#Dinh23a]_, based on [#Dinh21]_. Constant values are based on
    Table 3 of [#Dinh21]_ for proton exchange membrane (PEM) electrolysers
    predicted for the year 2030.

    .. math::
        AHP = \\frac{AEP}{\\frac{E_{electrolyser}}{\\eta} + E_{plant}}

    where :math:`AHP` is the annual hydrogen production [kg], :math:`AEP` is
    the annual energy production of wind farm [MWh], :math:`E_{electrolyser}`
    is the electricity required to supply the electrolyser to produce 1 kg of
    hydrogen [MWh kg⁻¹], :math:`\\eta` is the conversion efficiency of the
    electrolyser, and :math:`E_{plant}` is the electricity consumed by other
    parts of the hydrogen plant [MWh kg⁻¹].
    """
    return aep / (e_elec / eta_conv + e_pcl)


def electrolyser_capacity(n_turbines, wt_power=REF_RATED_POWER, loss=0.1):
    """
    Calculate the electrolyser capacity for a wind farm.

    Parameters
    ----------
    n_turbines : int
        Number of wind turbines in wind farm
    loss : float
        Electricity loss due to wake and transmission
    wt_power : float
        Rated power of the reference wind turbine [MW]

    Returns
    -------
    float
        Electrolyser capacity [MW]

    Notes
    -----
    Based on [#Dinh23b]_; the electrolyser capacity is the offshore wind farm
    capacity minus electricity loss. Electricity loss includes wake loss and
    transmission system loss. Since the electrolyser is situated at the
    offshore wind farm, only the wake loss is considered.
    """
    return n_turbines * wt_power * (1 - loss)


def capex_pipeline(e_cap, p_rate=0.0055, rho=8, u=15):
    """Capital expenditure (CAPEX) for the pipeline.

    Parameters
    ----------
    e_cap : float
        Electrolyser capacity [MW]
    p_rate : float
        Electrolyser production rate [kg s⁻¹ MW⁻¹]
    rho : float
        Mass density of hydrogen [kg m⁻³]
    u : float
        Average fluid velocity [m s⁻¹]

    Returns
    -------
    float
        CAPEX of the pipeline per km of pipeline [€ km⁻¹]

    Notes
    -----
    See Eqn. (18) of [#Dinh23a]_ and Section 3.1 of [#Dinh23c]_, from which
    the following text has been taken.

    The estimation of offshore pipeline costs is based on the onshore pipeline
    calculations of [#Baufume13]_, multiplied by a factor of two to reflect the
    estimated cost scaling of onshore to expected offshore costs, as suggested
    by [#IRENA22]_.

    In the reference case, the resulting capital expenditure for transmission
    pipelines and associated recompression stations was represented by a
    second order polynomial function.

    The electrolyser production rate was based on the Siemens - Silyzer 300
    [#Siemens]_.

    Because the electrolyser plant is assumed to be operating at full capacity
    at all times, the CAPEX was calculated considering a 75% utilisation rate,
    i.e. the pipeline capacity is 33% oversized [#IEA19]_.

    .. math::
        CAPEX = 2 \\, \\left( 16,000,000 \\, \\frac{P_{electrolyser}
        \\cdot EPR}{\\rho_{H2} \\cdot v_{H2} \\cdot \\pi} + 1,197,200 \\,
        \\sqrt{\\frac{P_{electrolyser} \\cdot EPR}{\\rho_{H2} \\cdot v_{H2}
        \\cdot \\pi}} + 329,000 \\right)

    where :math:`CAPEX` is the CAPEX of the pipeline per km of pipeline
    [€ km⁻¹], :math:`P_{electrolyser}` is the electrolyser capacity [MW],
    :math:`EPR` is the electrolyser production rate [kg s⁻¹ MW⁻¹],
    :math:`\\rho_{H2}` is the mass density of hydrogen [kg m⁻³], and
    :math:`v_{H2}` is the average fluid velocity [m s⁻¹].
    """
    f = e_cap * p_rate / (rho * u * np.pi)
    return 2e3 * (16000 * f + 1197.2 * np.sqrt(f) + 329)


def lcot_pipeline(
    capex,
    transmission_distance,
    ahp,
    opex_factor=0.02,
    discount_rate=0.08,
    lifetime=30,
):
    """Levelised cost of transmission (LCOT) of hydrogen in pipelines.

    Parameters
    ----------
    capex : float
        Capital expenditure (CAPEX) of the pipeline [€ km⁻¹]
    transmission_distance : float
        Pipeline transmission distance [km]
    ahp : float
        Annual hydrogen production [kg]
    opex_factor : float
        Ratio of the operational expenditure (OPEX) to the CAPEX; the OPEX is
        calculated as a percentage of the CAPEX
    discount_rate : float
        Discount rate
    lifetime : int
        Lifetime of the pipeline [year]

    Returns
    -------
    float
        LCOT of hydrogen in the pipeline [€ kg⁻¹]

    Notes
    -----
    See the introduction of Section 3, Eqn. (1) and (2), and Section 3.5, Eqn.
    (22) of [#Dinh23c]_; see Tables 2 and 3 for the assumptions and constants
    used.

    .. math::
        LCOT = \\frac{\\mathrm{total\\ lifecycle\\ costs\\ of\\ all\\
        components}}
        {\\mathrm{lifetime\\ hydrogen\\ transported}}
    .. math::
        LCOT = \\frac{CAPEX \\cdot d + \\sum_{l=0}^{L}
        \\frac{OPEX}{{(1 + DR)}^l}}
        {\\sum_{l=0}^{L} \\frac{AHP}{{(1 + DR)}^l}}

    where :math:`LCOT` is the LCOT of hydrogen in pipelines [€ kg⁻¹],
    :math:`CAPEX` is the CAPEX of the pipeline per km of pipeline
    [€ km⁻¹], :math:`d` is the pipeline transmission distance [km],
    :math:`OPEX` is the OPEX of the pipeline [€ kg⁻¹], :math:`DR` is the
    discount rate, :math:`AHP` is the annual hydrogen production [kg], and
    :math:`L` is the lifetime of the pipeline [year].
    """
    f = sum(
        1 / np.power((1 + discount_rate), year) for year in range(lifetime + 1)
    )
    opex = capex * opex_factor
    return (capex * transmission_distance + opex * f) / (ahp * f)


def rotor_area(diameter=REF_DIAMETER):
    """Wind turbine rotor swept area.

    Parameters
    ----------
    diameter : float
        Wind turbine rotor diameter [m]

    Returns
    -------
    float
        Wind turbine rotor swept area [m²]

    Notes
    -----
    .. math::
        A = \\frac{\\pi \\cdot D^2}{4}

    where :math:`A` is the wind turbine rotor swept area [m²] and :math:`D` is
    the rotor diameter [m].
    """
    return np.pi * np.square(diameter) / 4


def power_wind_resource(v, rho=1.225, diameter=REF_DIAMETER):
    """Total wind resource power passing through a wind turbine's rotor.

    Parameters
    ----------
    v : float
        Wind speed [m s⁻¹]
    rho : float
        Air density [kg m⁻³]
    diameter : float
        Wind turbine rotor diameter [m]

    Returns
    -------
    float
        Power contained in the wind resource [MW]

    Notes
    -----
    .. math::
        P_{wind} = \\frac{1}{2} \\,
        \\frac{\\rho_{air} \\cdot A \\cdot v^3}{10^6}

    where :math:`P_{wind}` is the power contained in the wind resource [MW],
    :math:`\\rho_{air}` is the air density [kg m⁻³], :math:`A` is the rotor
    swept area [m²], and :math:`v` is the wind speed [m s⁻¹].
    """
    return 0.5 * rho * rotor_area(diameter=diameter) * np.power(v, 3) / 1e6


def power_coefficient(v):
    """Power coefficient of the reference wind turbine.

    Parameters
    ----------
    v : float
        Wind speed [m s⁻¹]

    Returns
    -------
    float
        Power coefficient

    Notes
    -----
    This is specific to the power curve of the NREL 15 MW reference wind
    turbine [#Musial19]_.

    .. math::
        C_p = \\frac{P}{P_{wind}} = 2 \\, \\frac{P \\times 10^6}
        {\\rho_{air} \\cdot A \\cdot v^3}

    where :math:`P` is the wind turbine power output [MW], :math:`P_{wind}` is
    the power contained in the wind resource [MW]:math:`\\rho_{air}` is the
    air density [kg m⁻³], :math:`A` is the rotor swept area [m²], :math:`v` is
    the wind speed [m s⁻¹], and :math:`C_p` is the power coefficient of the
    wind turbine.
    """
    if REF_V_CUT_IN <= v <= REF_V_CUT_OUT:
        coeff = ref_power_curve(v=v) / power_wind_resource(v=v)
    else:
        coeff = 0
    return coeff
