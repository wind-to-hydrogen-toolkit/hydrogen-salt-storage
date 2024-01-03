"""test `src.optimisation` functions.

"""

import functools

import numpy as np
from scipy import integrate

from src import optimisation as opt


def test_ref_power_curve():
    """Test the `src.optimisation.ref_power_curve` function"""
    wind_speeds = list(range(31))
    power_curve = (
        [0] * 4
        + [0.499, 1.424, 2.732, 4.469, 6.643, 9.459, 12.975]
        + [15] * 15
        + [0] * 5
    )

    for v, p in zip(wind_speeds, power_curve):
        assert round(opt.ref_power_curve(v=v), 3) == p
        assert isinstance(opt.ref_power_curve(v=v), (float, int))


def test_weibull_probability_distribution():
    """Test the `src.optimisation.weibull_probability_distribution` function"""
    k_vals = [1.4, 1.7, 2.0, 2.15, 2.25]
    c_vals = [5.1, 9.2, 10.4, 8, 12]
    v_vals = [0, 7.2, 11.7, 14.6, 28, 20]
    weibull = []
    weibull_func = []
    for k, c, v in zip(k_vals, c_vals, v_vals):
        weibull.append(
            k / c * ((v / c) ** (k - 1)) * np.e ** (-((v / c) ** k))
        )
        weibull_func.append(
            opt.weibull_probability_distribution(k=k, c=c, v=v)
        )

    weibull = [round(x, 7) for x in weibull]
    weibull_func = [round(x, 7) for x in weibull_func]

    assert weibull_func == weibull


def test_annual_energy_production():
    """Test the `src.optimisation.annual_energy_production` function"""
    n_turbines = [50, 65, 80, 95]
    k_vals = [1.4, 1.7, 1.9, 2.0, 2.15]
    c_vals = [5.1, 9.2, 11, 10.4, 8]
    cut_in = 4
    cut_out = 25
    w_loss = 0.1
    aep = []
    aep_func = []
    for n, k, c in zip(n_turbines, k_vals, c_vals):
        aep_func.append(
            opt.annual_energy_production(n_turbines=n, k=k, c=c)[0]
        )
        integration = integrate.quad(
            functools.partial(opt.weibull_power_curve, k=k, c=c),
            a=cut_in,
            b=cut_out,
            limit=100,  # an upper bound on the number of subintervals used
            epsabs=1.49e-6,  # absolute error tolerance
        )
        aep.append(365 * 24 * n * (1 - w_loss) * integration[0])

    assert aep_func == aep


def test_lcot_pipeline():
    """Test the `src.optimisation.lcot_pipeline` function"""
    capex = 1000
    transmission_distance = 100
    ahp = 500
    opex_factor = 0.03
    discount_rate = 0.05
    lifetime = 40
    opex = capex * opex_factor

    lcot = (
        capex * transmission_distance
        + sum(
            opex / np.power((1 + discount_rate), year)
            for year in range(lifetime + 1)
        )
    ) / sum(
        ahp / np.power((1 + discount_rate), year)
        for year in range(lifetime + 1)
    )
    lcot_func = opt.lcot_pipeline(
        capex=capex,
        transmission_distance=transmission_distance,
        ahp=ahp,
        opex_factor=opex_factor,
        discount_rate=discount_rate,
        lifetime=lifetime,
    )

    assert lcot_func == lcot
