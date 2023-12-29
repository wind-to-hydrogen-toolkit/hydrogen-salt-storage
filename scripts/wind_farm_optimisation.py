#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src import optimisation as opt

opt.rotor_area()

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

# opt.annual_energy_production(n_turbines=25, k=2, c=10.9)

fig, axes = plt.subplots(1, 1, figsize=(12, 6))
sns.lineplot(
    data=ref_data,
    x="wind_speed",
    y="power_curve",
    legend=False,
    color="crimson",
)
axes.grid(axis="y")
axes.set_xlim(0, 30)
sns.despine()
axes.set_xlabel("Wind speed [m s-1]")
axes.set_ylabel("Power [MW]")
plt.show()

fig, axes = plt.subplots(1, 1, figsize=(12, 6))
sns.lineplot(
    data=ref_data, x="wind_speed", y="weibull", legend=False, color="crimson"
)
axes.grid(axis="y")
axes.set_xlim(0, 30)
sns.despine()
axes.set_xlabel("Wind speed [m s-1]")
axes.set_ylabel("Weibull probability distribution function [s m-1]")
plt.show()
