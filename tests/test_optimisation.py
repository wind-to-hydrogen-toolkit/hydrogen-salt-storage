"""test `src.optimisation` functions.

"""

import numpy as np

from src import optimisation as opt


def test_ref_power_curve():
    """Test the `src.optimisation.ref_power_curve` function."""
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


def test_lcot_pipeline():
    """Test the `src.optimisation.lcot_pipeline` function."""
    capex = 1000
    transmission_distance = 100
    prod_h2 = 500
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
        prod_h2 / np.power((1 + discount_rate), year)
        for year in range(lifetime + 1)
    )
    lcot_func = opt.lcot_pipeline(
        capex=capex,
        transmission_distance=transmission_distance,
        prod_h2=prod_h2,
        opex_factor=opex_factor,
        discount_rate=discount_rate,
        lifetime=lifetime,
    )

    assert lcot_func == lcot
