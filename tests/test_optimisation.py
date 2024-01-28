"""test ``h2ss.optimisation`` functions.

"""

import numpy as np
from scipy import integrate

from h2ss import optimisation as opt


def test_ref_power_curve():
    """Test ``h2ss.optimisation.ref_power_curve``"""
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
        power_coeff = opt.power_coefficient(v=v) * opt.power_wind_resource(v=v)
        assert round(power_coeff, 10) == round(opt.ref_power_curve(v=v), 10)


def test_weibull_probability_distribution():
    """Test ``h2ss.optimisation.weibull_probability_distribution``"""
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
    weibull = [round(x, 10) for x in weibull]
    weibull_func = [round(x, 10) for x in weibull_func]
    assert weibull_func == weibull


def test_number_of_turbines():
    """Test ``h2ss.optimisation.number_of_turbines``"""
    owf_cap_list = [500, 600, 700, 800, 900]
    wt_power_list = [10, 11, 12, 13, 14]
    n_turbines = []
    n_turbines_func = []
    for owf_cap, wt_power in zip(owf_cap_list, wt_power_list):
        n_turbines_func.append(
            opt.number_of_turbines(owf_cap=np.array(owf_cap),
            wt_power=np.array(wt_power))
        )
        n_turbines.append(int(owf_cap / wt_power))
    assert n_turbines_func == n_turbines


def integrate_lambda(k, c):
    """Lambda function for the quad integration test below"""
    return lambda v: (
        opt.ref_power_curve(v=v)
        * opt.weibull_probability_distribution(k=k, c=c, v=v)
    )


def test_annual_energy_production():
    """Test ``h2ss.optimisation.annual_energy_production``"""
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
            integrate_lambda(k=k, c=c),
            a=cut_in,
            b=cut_out,
            limit=100,  # an upper bound on the number of subintervals used
            epsabs=1.49e-6,  # absolute error tolerance
        )
        aep.append(365 * 24 * n * (1 - w_loss) * integration[0])
    assert aep_func == aep


def test_annual_hydrogen_production():
    """Test ``h2ss.optimisation.annual_hydrogen_production``"""
    aep = 15e6
    e_elec = 0.07
    eta_conv = 0.9
    e_pcl = 0.004
    assert opt.annual_hydrogen_production(
        aep=aep, e_elec=e_elec, eta_conv=eta_conv, e_pcl=e_pcl
    ) == aep / (e_elec / eta_conv + e_pcl)


def test_electrolyser_capacity():
    """Test ``h2ss.optimisation.electrolyser_capacity``"""
    n_turbines_list = [25, 50, 75, 100, 125]
    cap_ratio_list = [0.9, 0.88, 0.86, 0.84, 0.82]
    e_cap = []
    e_cap_func = []
    for n_turbines, cap in zip(n_turbines_list, cap_ratio_list):
        e_cap_func.append(
            opt.electrolyser_capacity(n_turbines=n_turbines, cap_ratio=cap)
        )
        owf_cap = n_turbines * opt.REF_RATED_POWER
        e_cap.append(owf_cap * cap)
    assert e_cap_func == e_cap


def test_capex_pipeline():
    """Test ``h2ss.optimisation.capex_pipeline``"""
    e_cap = 500
    p_rate = 0.006
    capex = 2 * (
        16000000 * e_cap * p_rate / (8 * 15 * np.pi)
        + 1197200 * np.sqrt(e_cap * p_rate / (8 * 15 * np.pi))
        + 329000
    )
    assert opt.capex_pipeline(e_cap=e_cap, p_rate=p_rate) == capex


def test_lcot_pipeline():
    """Test ``h2ss.optimisation.lcot_pipeline``"""
    capex = 1000
    transmission_distance = 100
    ahp = 500
    opex_ratio = 0.03
    discount_rate = 0.05
    lifetime = 40
    opex = capex * opex_ratio
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
        opex_ratio=opex_ratio,
        discount_rate=discount_rate,
        lifetime=lifetime,
    )
    assert lcot_func == lcot
