"""Test `src.capacity` functions.

.. rubric:: References
.. [#Caglayan20] Caglayan, D. G., Weber, N., Heinrichs, H. U., Linßen, J.,
    Robinius, M., Kukla, P. A., and Stolten, D. (2020). ‘Technical potential
    of salt caverns for hydrogen storage in Europe’, International Journal of
    Hydrogen Energy, 45(11), pp. 6793–6805.
    https://doi.org/10.1016/j.ijhydene.2019.12.161.
"""

import numpy as np
import pandas as pd
from pandas.testing import assert_series_equal
from pyfluids import Fluid, FluidsList, Input

from src import capacity as cap


def test_cavern_volume():
    """Test the `src.capacity.cavern_volume` function"""
    h_cylinder = [10, 20, 30, 40, 50]
    h_cone = [2, 4, 6, 8, 10]
    roof_angle = [5, 10, 15, 20, 25]
    v_cavern = []
    v_cavern_func = []

    for cyl, con, a in zip(h_cylinder, h_cone, roof_angle):
        r = con / np.tan(np.deg2rad(a))
        v_cylinder = np.pi * np.square(r) * cyl
        v_cone = np.pi * np.square(r) * con / 3
        v_cavern.append(v_cylinder + v_cone * 2)
        v_cavern_func.append(
            cap.cavern_volume(height=cyl + con * 2, diameter=r * 2, theta=a)
        )

    v_cavern = [round(x, 5) for x in v_cavern]
    v_cavern_func = [round(x, 5) for x in v_cavern_func]

    assert v_cavern_func == v_cavern


def test_corrected_cavern_volume():
    """Test the `src.capacity.corrected_cavern_volume` function"""
    correction_factors = 0.7 * (1 - 0.25 * 0.865 * 1.46)
    volumes = [1e4, 2.5e4, 3e4, 4.7e4, 5.3e4]
    v_cavern = []
    v_cavern_func = []

    for v in volumes:
        v_cavern.append(v * correction_factors)
        v_cavern_func.append(cap.corrected_cavern_volume(v_cavern=v))

    assert v_cavern_func == v_cavern


def test_temperature_cavern_mid_point():
    """Test the `src.capacity.temperature_cavern_mid_point` function"""
    heights = [50, 100, 150, 200, 250]
    top_depths = [600, 800, 1100, 1300, 1700]
    t_delta = [35, 37, 39, 41, 43]
    t_0 = [8, 9, 10, 11, 12]
    t_mid_point = []
    t_mid_point_func = []

    for h, d, td, t0 in zip(heights, top_depths, t_delta, t_0):
        t_mid_point.append(t0 + td * (d + 0.5 * h) / 1000 + 273.15)
        t_mid_point_func.append(
            cap.temperature_cavern_mid_point(
                height=h, depth_top=d, t_0=t0, t_delta=td
            )
        )

    assert t_mid_point_func == t_mid_point


def test_pressure_operating():
    """Test the `src.capacity.pressure_operating` function"""
    thickness_overburden = [550, 650, 750, 850, 950]
    p_operating = []
    p_operating_func = []

    for t in thickness_overburden:
        p_casing = (2400 * t + 2200 * 50) * 9.81
        p_operating.append((0.3 * p_casing, 0.8 * p_casing))
        p_operating_func.append(cap.pressure_operating(thickness_overburden=t))

    assert p_operating_func == p_operating


def test_density_hydrogen_gas():
    """Test the `src.capacity.density_hydrogen_gas` function.

    Notes
    -----
    Use Eqn. (3) of [#Caglayan20]_ to derive the density and compare it with
    the density obtained from the `src.capacity.density_hydrogen_gas` function.
    The values should be approximately the same (rounded to one decimal place).
    PyFluids is used to just derive the compressibility factor, :math:`Z`, for
    the former. The `pyfluids.ini` configuration file has been set such that
    the default units are SI units.

    .. math::
        \\rho = \\frac{P \\cdot M}{Z \\cdot R \\cdot T}

    where :math:`M` is the molar mass of hydrogen gas [kg mol⁻¹], :math:`R` is
    the universal gas constant [J K⁻¹ mol⁻¹], :math:`P` is the pressure [Pa],
    :math:`T` is the temperature [K], and :math:`\\rho` is the hydrogen gas
    density [kg m⁻³].
    """
    p_operating_min = [2.03e7, 3.29e6, 3.29e6, 5.5e6, 1.5e6]
    p_operating_max = [3.698e7, 2.98e7, 7.62e6, 3.9e7, 2e7]
    t_mid_point = [327.4, 303.1, 362.75, 300, 358.9]

    m = 0.00201588  # molar mass of hydrogen gas [kg mol⁻¹]
    r = 8.314  # universal gas constant [J K⁻¹ mol⁻¹]

    rho_h2 = []

    for p_min, p_max, t in zip(p_operating_min, p_operating_max, t_mid_point):
        h2_min = Fluid(FluidsList.Hydrogen).with_state(
            Input.pressure(p_min), Input.temperature(t)
        )
        h2_max = Fluid(FluidsList.Hydrogen).with_state(
            Input.pressure(p_max), Input.temperature(t)
        )

        rho_approx_min = p_min * m / (h2_min.compressibility * r * t)
        rho_approx_max = p_max * m / (h2_max.compressibility * r * t)

        rho_h2.append((rho_approx_min, rho_approx_max))

    rho_h2_func = cap.density_hydrogen_gas(
        p_operating_min=p_operating_min,
        p_operating_max=p_operating_max,
        t_mid_point=t_mid_point,
    )

    rho_h2 = pd.DataFrame(rho_h2)

    assert_series_equal(round(rho_h2_func[0], 1), round(rho_h2[0], 1))
    assert_series_equal(round(rho_h2_func[1], 1), round(rho_h2[1], 1))


def test_mass_hydrogen_working():
    """Test the `src.capacity.mass_hydrogen_working` function"""
    rho_h2_min = [2.7, 8.5, 6, 4.563, 3.5]
    rho_h2_max = [20, 12, 16, 6.6, 10.4928]
    v_cavern = [2.7e5, 3.5e5, 4.1234e5, 6e5, 7.2e5]
    m_working = []
    m_working_func = []

    for mi, ma, v in zip(rho_h2_min, rho_h2_max, v_cavern):
        m_working.append((ma - mi) * v)
        m_working_func.append(
            cap.mass_hydrogen_working(rho_h2_min=mi, rho_h2_max=ma, v_cavern=v)
        )

    assert m_working_func == m_working


def test_energy_storage_capacity():
    """Test the `src.capacity.energy_storage_capacity` function"""
    m_working = [6.4e5, 1.5e6, 2.6e6, 8.2e6, 3.5e6]
    capacity = []
    capacity_func = []

    for m in m_working:
        capacity.append(m * 119.96 / 3600000)
        capacity_func.append(cap.energy_storage_capacity(m_working=m))

    assert capacity_func == capacity
