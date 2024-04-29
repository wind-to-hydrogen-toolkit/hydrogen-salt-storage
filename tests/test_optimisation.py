"""test ``h2ss.optimisation`` functions.

"""

import numpy as np
from scipy import integrate
import pandas as pd
from pandas.testing import assert_frame_equal
from geopandas.testing import assert_geodataframe_equal, assert_geoseries_equal
import geopandas as gpd
from shapely.geometry import Point, Polygon

from h2ss import optimisation as opt
from h2ss import data as rd


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
            opt.number_of_turbines(
                owf_cap=np.array(owf_cap), wt_power=np.array(wt_power)
            )
        )
        n_turbines.append(int(owf_cap / wt_power))
    assert n_turbines_func == n_turbines


def integrate_lambda(k, c):
    """Lambda function for the quad integration test below"""
    return lambda v: (
        opt.ref_power_curve(v=v)
        * opt.weibull_probability_distribution(k=k, c=c, v=v)
    )


def test_annual_energy_production_function():
    """Test ``h2ss.optimisation.annual_energy_production_function``"""
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
            opt.annual_energy_production_function(n_turbines=n, k=k, c=c)[0]
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


def test_transmission_distance():
    """Test ``h2ss.optimisation.transmission_distance``"""
    wf_data = gpd.GeoDataFrame({"name": ["Dublin Array", "Codling", "NISA"], "geometry": [Polygon([(104.0, -27.0), (104.0, -26.9), (104.1, -26.9), (104.1, -27.0), (104.0, -27.0)]), Polygon([(105.0, -26.0), (105.0, -25.9), (105.1, -25.9), (105.1, -26.0), (105.0, -26.0)]), Polygon([(103.0, -28.0), (103.0, -27.9), (103.1, -27.9), (103.1, -28.0), (103.0, -28.0)])]}, crs=4326).to_crs(rd.CRS)
    cavern_df = gpd.GeoDataFrame({"geometry": [Point([111.0, -29.0]), Point([111.1, -29.1]), Point([111.2, -29.2]), Point([111.3, -29.3]), Point([111.4, -29.4])]}, crs=4326).to_crs(rd.CRS)
    injection_point_coords = (114, 10, -28, -1)
    lond, lonm, latd, latm = injection_point_coords
    injection_point = gpd.GeoSeries(
        [Point(((lond) + (lonm) / 60), ((latd) + (latm) / 60))], crs=4326
    ).to_crs(rd.CRS)
    distance_ip = []
    for j in range(len(cavern_df)):
        distance_ip.append(
            injection_point.distance(
                cavern_df.iloc[[j]]["geometry"], align=False
            ).values[0]
            / 1000
        )
    cavern_df["distance_ip"] = distance_ip
    distance_wf = {}
    for i in range(len(wf_data)):
        distance_wf[wf_data["name"][i]] = []
        for j in range(len(cavern_df)):
            distance_wf[wf_data["name"][i]].append(
                (
                    wf_data.iloc[[i]]
                    .distance(cavern_df.iloc[[j]]["geometry"], align=False)
                    .values[0]
                    / 1000
                    + cavern_df.iloc[[j]]["distance_ip"].values[0]
                )
            )
        cavern_df[f"dist_{wf_data['name'][i].replace(' ', '_')}"] = (
            distance_wf[wf_data["name"][i]]
        )
    cavern_df_func, injection_point_func = opt.transmission_distance(cavern_df=cavern_df, wf_data=wf_data, injection_point_coords=injection_point_coords)
    assert_geodataframe_equal(cavern_df_func, cavern_df)
    assert_geoseries_equal(injection_point_func, injection_point)


def test_electrolyser_capacity():
    """Test ``h2ss.optimisation.electrolyser_capacity``"""
    n_turbines_list = [25, 50, 75, 100, 125]
    cap_ratio_list = [0.9, 0.88, 0.86, 0.84, 0.82]
    e_cap = []
    e_cap_func = []
    for n_turbines, cap in zip(n_turbines_list, cap_ratio_list):
        e_cap_func.append(
            opt.electrolyser_capacity(
                n_turbines=np.array(n_turbines), cap_ratio=np.array(cap)
            )
        )
        owf_cap = n_turbines * opt.REF_RATED_POWER
        e_cap.append(int(owf_cap * cap))
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


def test_lcot_pipeline_function():
    """Test ``h2ss.optimisation.lcot_pipeline_function``"""
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
    lcot_func = opt.lcot_pipeline_function(
        capex=capex,
        d_transmission=transmission_distance,
        ahp=ahp,
        opex_ratio=opex_ratio,
        discount_rate=discount_rate,
        lifetime=lifetime,
    )
    assert lcot_func == lcot


def test_lcot_pipeline():
    """Test ``h2ss.optimisation.lcot_pipeline``"""
    weibull_wf_data = pd.DataFrame({"name": ["Dublin Array", "Codling", "NISA"], "CAPEX": [1000, 2000, 3000], "AHP": [8000, 9000, 10000]})
    cavern_df = pd.DataFrame({"dist_Dublin_Array": [x + 2 for x in range(20)], "dist_Codling": [x + 3 for x in range(20)], "dist_NISA": [x + 4 for x in range(20)]})

    for wf, capex, ahp in zip(list(weibull_wf_data["name"]), list(weibull_wf_data["CAPEX"]), list(weibull_wf_data["AHP"])):
        cavern_df[f"LCOT_{wf.replace(' ', '_')}"] = opt.lcot_pipeline_function(
            capex=capex,
            d_transmission=cavern_df[f"dist_{wf.replace(' ', '_')}"],
            ahp=ahp,
        )

    cavern_df_func = opt.lcot_pipeline(weibull_wf_data=weibull_wf_data, cavern_df=cavern_df)
    assert_frame_equal(cavern_df_func, cavern_df)
