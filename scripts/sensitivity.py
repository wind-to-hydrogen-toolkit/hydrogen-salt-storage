#!/usr/bin/env python
# coding: utf-8

# # Sensitivity analysis

import glob
import os
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from h2ss import compare

cavern_diameter = np.arange(80, 101, step=1)
cavern_height = np.arange(85, 312, step=1)


def generate_sensitivity_data(cavern_diameter, cavern_height):
    """Generate data to perform sensitivity analysis"""
    os.makedirs("data", "sensitivity", exist_ok=True)
    ds, extent, exclusions = compare.load_all_data()
    for d, h in product(cavern_diameter, cavern_height):
        df = compare.capacity_function(ds, extent, exclusions, d, h)
        df.to_csv(
            os.path.join("data", "sensitivity", f"sensitivity_d{d}_h{h}.csv")
        )
        print(f"sensitivity_d{d}_h{h}.csv done!")


# generate_sensitivity_data(cavern_diameter, cavern_height)

len(list(product(cavern_diameter, cavern_height)))

len(glob.glob(os.path.join("data", "sensitivity", f"sensitivity_*.csv")))

df = pd.concat(
    (
        pd.read_csv(f)
        for f in glob.glob(
            os.path.join("data", "sensitivity", f"sensitivity_*.csv")
        )
    ),
    ignore_index=True,
)

df.drop(columns=["Unnamed: 0"], inplace=True)

df.describe()

sns.scatterplot(
    data=df,
    hue="cavern_diameter",
    y="capacity",
    x="cavern_height",
    linewidth=0,
    alpha=0.75,
    s=5,
)
sns.despine()
plt.show()

# ## Cavern height - mean capacity

data = df.groupby(["cavern_height"]).sum()[["capacity"]].reset_index()
sns.lineplot(data=data, x="cavern_height", y="capacity")
sns.despine()
plt.show()

# ## Cavern height - total capacity

data = df.groupby(["cavern_height"]).mean()[["capacity"]].reset_index()
sns.lineplot(data=data, x="cavern_height", y="capacity")
sns.despine()
plt.show()

# ## Cavern diameter - mean capacity

data = df.groupby(["cavern_diameter"]).sum()[["capacity"]].reset_index()
sns.lineplot(data=data, x="cavern_diameter", y="capacity")
sns.despine()
plt.show()

# ## Cavern diameter - total capacity

data = df.groupby(["cavern_diameter"]).mean()[["capacity"]].reset_index()
sns.lineplot(data=data, x="cavern_diameter", y="capacity")
sns.despine()
plt.show()

# ## Base case

base = df[
    (df["cavern_diameter"] == 85) & (df["cavern_height"] == 120)
].reset_index(drop=True)

base.describe()[["capacity"]]

base_mean = base[["capacity"]].mean().values[0]

base_sum = base[["capacity"]].sum().values[0]

print(f"{base_sum:.4f}")

# ## Base diameter, varying height, mean capacity

dd = df[(df["cavern_diameter"] == 85)].reset_index(drop=True)

dd_diff = (
    pd.DataFrame(dd.groupby("cavern_height").mean()["capacity"] - base_mean)
    / base_mean
    * 100
)

plt.figure(figsize=(8, 5))
ax = sns.scatterplot(
    data=dd_diff,
    hue="cavern_height",
    y="capacity",
    x="cavern_height",
    palette="icefire_r",
    legend=False,
    linewidth=0,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline(120, color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern height [m]")
ax.set_ylabel("Difference in mean\nhydrogen storage capacity [%]")
sns.despine()
plt.show()

# ## Base diameter, varying height, total capacity

dd_diff = (
    pd.DataFrame(dd.groupby("cavern_height").sum()["capacity"] - base_sum)
    / base_sum
    * 100
)

plt.figure(figsize=(8, 5))
ax = sns.scatterplot(
    data=dd_diff,
    hue="cavern_height",
    y="capacity",
    x="cavern_height",
    palette="icefire",
    legend=False,
    linewidth=0,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline(120, color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern height [m]")
ax.set_ylabel("Difference in total\nhydrogen storage capacity [%]")
sns.despine()
plt.show()

# ## Base height, varying diameter, mean capacity

dh = df[(df["cavern_height"] == 120)].reset_index(drop=True)

dh_diff = (
    pd.DataFrame(dh.groupby("cavern_diameter").mean()["capacity"] - base_mean)
    / base_mean
    * 100
)

plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=dh_diff,
    hue="cavern_diameter",
    y="capacity",
    x="cavern_diameter",
    palette="icefire_r",
    legend=False,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline("85", color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Difference in mean\nhydrogen storage capacity [%]")
sns.despine()
plt.show()

# ## Base height, varying diameter, total capacity

dh_diff = (
    pd.DataFrame(dh.groupby("cavern_diameter").sum()["capacity"] - base_sum)
    / base_sum
    * 100
)

plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=dh_diff,
    hue="cavern_diameter",
    y="capacity",
    x="cavern_diameter",
    palette="icefire",
    legend=False,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline("85", color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Difference in total\nhydrogen storage capacity [%]")
sns.despine()
plt.show()
