"""test_functions.py

.. rubric:: References
.. [#Caglayan20] Caglayan, D. G., Weber, N., Heinrichs, H. U., Linßen, J.,
    Robinius, M., Kukla, P. A., and Stolten, D. (2020). ‘Technical potential
    of salt caverns for hydrogen storage in Europe’, International Journal of
    Hydrogen Energy, 45(11), pp. 6793–6805.
    https://doi.org/10.1016/j.ijhydene.2019.12.161.
"""

# import os

import pandas as pd
from pandas.testing import assert_series_equal
from pyfluids import Fluid, FluidsList, Input

from src import capacity as cap


def test_density_hydrogen_gas():
    """Test the `capacity.density_hydrogen_gas` function.

    Notes
    -----
    Use Eqn. (3) of [#Caglayan20]_ to derive the density and compare it with
    the density obtained from the `capacity.density_hydrogen_gas` function.
    PyFluids is used to just derive the compressibility factor, :math:`Z`, for
    the former. :math:`M` is the molar mass of hydrogen gas [kg mol⁻¹],
    :math:`R` is the universal gas constant [J K⁻¹ mol⁻¹], :math:`P` is the
    pressure [Pa], :math:`T` is the temperature [K], and :math:`\\rho` is the
    hydrogen gas density [kg m⁻³]. The `pyfluids.ini` configuration file has
    been set such that the default units are SI units.

    .. math::
        \\rho = \\frac{P \\times M}{Z \\times R \\times T}
    """
    p_operating_min = [2.03e7, 3.698e7, 3.29e6, 3e7, 5.5e6]
    p_operating_max = [3.698e7, 3.29e6, 7.62e6, 5.5e6, 2e7]
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
