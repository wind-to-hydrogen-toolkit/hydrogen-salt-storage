#!/usr/bin/env python
# coding: utf-8

# # Wind farm optimisation

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src import optimisation as opt

ref_data = {}

ref_data["wind_speed"] = [0 + 0.01 * n for n in range(3100)]
ref_data["power_curve"] = []
ref_data["weibull"] = []

for v in ref_data["wind_speed"]:
    ref_data["power_curve"].append(opt.ref_power_curve(v=v))
    ref_data["weibull"].append(
        opt.weibull_probability_distribution(v=v, k=2, c=10.9)
    )

ref_data = pd.DataFrame(ref_data)

ref_data.head()

opt.annual_energy_production(n_turbines=40, k=2, c=10.9)

ax = ref_data.plot(
    x="wind_speed",
    y="power_curve",
    label="Power",
    ylabel="Power [MW]",
    linewidth=2,
    color="crimson",
    figsize=(12, 6),
)
ref_data.plot(
    x="wind_speed",
    y="weibull",
    label="Weibull",
    ylabel=(
        "Weibull probability distribution function "
        + "[s m\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]"
    ),
    linewidth=2,
    color="royalblue",
    linestyle="dashed",
    secondary_y=True,
    ax=ax,
)
ax.set_xlabel("Wind speed [m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")
sns.despine(top=True, right=False)
