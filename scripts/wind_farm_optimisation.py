#!/usr/bin/env python
# coding: utf-8

# # Wind farm optimisation

import importlib
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src import functions as fns
from src import optimisation as opt

importlib.reload(opt)

# ## Power curve [MW] and Weibull wind speed distribution

ds, extent = fns.read_dat_file(dat_path=os.path.join("data", "kish-basin"))

# extract data for wind farms at 150 m
weibull = fns.read_weibull_data(
    data_path_weibull=os.path.join(
        "data", "weibull-parameters-wind-speeds", "Weibull_150m_params_ITM.zip"
    ),
    data_path_wind_farms=os.path.join(
        "data",
        "wind-farms-foreshore-process",
        "wind-farms-foreshore-process.zip",
    ),
    dat_extent=extent,
)

weibull

ref_data = {}

# generate Weibull distribution
for n in weibull["Name"]:
    ref_data[n] = {}
    ref_data[n]["wind_speed"] = [0 + 0.01 * n for n in range(3100)]
    ref_data[n]["power_curve"] = []
    ref_data[n][n] = []
    for v in ref_data[n]["wind_speed"]:
        ref_data[n]["power_curve"].append(opt.ref_power_curve(v=v))
        ref_data[n][n].append(
            opt.weibull_probability_distribution(
                v=v,
                k=weibull[weibull["Name"] == n][("k", "mean")].iloc[0],
                c=weibull[weibull["Name"] == n][("c", "mean")].iloc[0],
            )
        )
    ref_data[n] = pd.DataFrame(ref_data[n])

ref_data = pd.concat(ref_data.values(), axis=1).T.drop_duplicates().T

ref_data.head()

ax = ref_data.plot(
    x="wind_speed",
    y="power_curve",
    ylabel="Power [MW]",
    linewidth=3,
    color="crimson",
    figsize=(12, 6),
    legend=False,
)
ax.set_xlabel("Wind speed [m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")
ax.set_ylabel("Power [MW]")
sns.despine()

ax = ref_data.drop(columns=["power_curve"]).plot(
    x="wind_speed",
    cmap="flare_r",
    figsize=(12, 6),
    linestyle="dashed",
    linewidth=2,
    alpha=0.85,
)
ax.set_xlabel("Wind speed [m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]")
ax.set_ylabel(
    "Weibull probability distribution function "
    "[s m\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]"
)
sns.despine()

# ## Annual energy production [MWh]

# max wind farm capacity
weibull["capacity"] = [1300, 824, 500]

# number of 15 MW turbines, rounded down to the nearest integer
weibull["n_turbines"] = (weibull["capacity"] / 15).astype(int)

aep = []
integral = []
abserr = []
for n in weibull["Name"]:
    aepwt = opt.annual_energy_production(
        n_turbines=weibull[weibull["Name"] == n]["n_turbines"].iloc[0],
        k=weibull[weibull["Name"] == n][("k", "mean")].iloc[0],
        c=weibull[weibull["Name"] == n][("c", "mean")].iloc[0],
    )
    aep.append(aepwt)

aep = pd.DataFrame(aep)
aep.columns = ["AEP", "integral", "abserr"]

aep = pd.concat([weibull, aep], axis=1)

aep

# ## Annual hydrogen production [kg]

aep["AHP"] = opt.annual_hydrogen_production(aep=aep["AEP"])

aep
