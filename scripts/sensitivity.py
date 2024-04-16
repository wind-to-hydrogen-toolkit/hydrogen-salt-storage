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

df["cavern_height"] = df["cavern_height"].astype(int)

df.describe()

ax = sns.scatterplot(
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

# ## Mean capacity

data = (
    df.groupby(["cavern_height", "cavern_diameter"])
    .mean()
    .reset_index()
    .pivot(index="cavern_height", columns="cavern_diameter", values="capacity")
    .sort_index(ascending=False)
)
f, ax = plt.subplots(figsize=(9, 7))
sns.heatmap(
    data,
    ax=ax,
    cmap="rocket_r",
    cbar_kws={"label": "Mean capacity [GWh]"},
)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Cavern height [m]")
ax.tick_params(axis="y", labelsize=11)
ax.tick_params(axis="x", labelsize=11)
plt.show()

# ## Total capacity

data = df.copy()
data["capacity"] = data["capacity"] / 1000
data = (
    data.groupby(["cavern_height", "cavern_diameter"])
    .sum()
    .reset_index()
    .pivot(index="cavern_height", columns="cavern_diameter", values="capacity")
    .sort_index(ascending=False)
)
f, ax = plt.subplots(figsize=(9, 7))
sns.heatmap(
    data, ax=ax, cmap="rocket_r", cbar_kws={"label": "Total capacity [TWh]"}
)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Cavern height [m]")
ax.tick_params(axis="y", labelsize=11)
ax.tick_params(axis="x", labelsize=11)
plt.show()

# ## Base case

base = df[
    (df["cavern_diameter"] == 85) & (df["cavern_height"] == 120)
].reset_index(drop=True)

base.describe()[["capacity"]]

base_mean = base[["capacity"]].mean().values[0]

base_sum = base[["capacity"]].sum().values[0]

print(f"{base_sum:.3f}")

# ## Base diameter, varying height

dd = df[(df["cavern_diameter"] == 85)].reset_index(drop=True)

dd_mean = (
    pd.DataFrame(dd.groupby("cavern_height").mean()["capacity"] - base_mean)
    / base_mean
    * 100
).reset_index()

plt.figure(figsize=(8, 5))
ax = sns.scatterplot(
    data=dd_mean,
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
ax.set_ylabel("Difference in mean storage capacity [%]")
sns.despine()
plt.show()

dd_mean1 = dd_mean[
    dd_mean["cavern_height"].isin([90 + 10 * n for n in range(21)])
]

plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=dd_mean1,
    hue="cavern_height",
    y="capacity",
    x="cavern_height",
    palette="icefire_r",
    legend=False,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline("120", color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern height [m]")
ax.set_ylabel("Difference in mean storage capacity [%]")
sns.despine()
plt.show()

dd_sum = (
    pd.DataFrame(dd.groupby("cavern_height").sum()["capacity"] - base_sum)
    / base_sum
    * 100
).reset_index()

plt.figure(figsize=(8, 5))
ax = sns.scatterplot(
    data=dd_sum,
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
ax.set_ylabel("Difference in total storage capacity [%]")
sns.despine()
plt.show()

dd_sum1 = dd_sum[
    dd_sum["cavern_height"].isin([90 + 10 * n for n in range(21)])
]

plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=dd_sum1,
    hue="cavern_height",
    y="capacity",
    x="cavern_height",
    palette="icefire",
    legend=False,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline("120", color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern height [m]")
ax.set_ylabel("Difference in total storage capacity [%]")
sns.despine()
plt.show()

# ## Base height, varying diameter

dh = df[(df["cavern_height"] == 120)].reset_index(drop=True)

dh_mean = (
    pd.DataFrame(dh.groupby("cavern_diameter").mean()["capacity"] - base_mean)
    / base_mean
    * 100
).reset_index()

plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=dh_mean,
    hue="cavern_diameter",
    y="capacity",
    x="cavern_diameter",
    palette="icefire_r",
    legend=False,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline("85", color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Difference in mean storage capacity [%]")
sns.despine()
plt.show()

dh_sum = (
    pd.DataFrame(dh.groupby("cavern_diameter").sum()["capacity"] - base_sum)
    / base_sum
    * 100
).reset_index()

plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=dh_sum,
    hue="cavern_diameter",
    y="capacity",
    x="cavern_diameter",
    palette="icefire",
    legend=False,
)
ax.axhline(0, color="darkslategrey", linewidth=1)
ax.axvline("85", color="darkslategrey", linewidth=1)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Difference in total storage capacity [%]")
sns.despine()
plt.show()

# ## Combined plots

dh_sum["type"] = "Total for the Kish Basin"
dh_mean["type"] = "Mean for a single cavern"
dd_sum["type"] = "Total for the Kish Basin"
dd_mean["type"] = "Mean for a single cavern"
dd_sum1["type"] = "Total for the Kish Basin"
dd_mean1["type"] = "Mean for a single cavern"

f, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
sns.scatterplot(
    data=pd.concat([dh_sum, dh_mean]),
    y="capacity",
    x="cavern_diameter",
    hue="type",
    linewidth=0,
    alpha=0.75,
    palette=sns.color_palette(["tab:blue", "tab:red"]),
    ax=ax[0],
)
sns.scatterplot(
    data=pd.concat([dd_sum, dd_mean]),
    y="capacity",
    x="cavern_height",
    hue="type",
    linewidth=0,
    alpha=0.75,
    palette=sns.color_palette(["tab:blue", "tab:red"]),
    ax=ax[1],
    legend=False,
)
for a in ax.flat:
    a.axhline(0, color="darkslategrey", linewidth=1)
    a.xaxis.grid(True, linewidth=0.25)
    a.yaxis.grid(True, linewidth=0.25)
ax[0].axvline(85, color="darkslategrey", linewidth=1)
ax[1].axvline(120, color="darkslategrey", linewidth=1)
ax[0].set_xlabel("Cavern diameter [m]")
ax[1].set_xlabel("Cavern height [m]")
ax[0].set_ylabel("Difference in storage capacity [%]")
ax[0].legend(title=None)
sns.despine()
plt.tight_layout()
plt.show()

f, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
sns.barplot(
    data=pd.concat([dh_sum, dh_mean]),
    hue="type",
    y="capacity",
    x="cavern_diameter",
    palette=sns.color_palette(["tab:blue", "tab:red"]),
    ax=ax[0],
)
sns.barplot(
    data=pd.concat([dd_sum1, dd_mean1]),
    hue="type",
    y="capacity",
    x="cavern_height",
    palette=sns.color_palette(["tab:blue", "tab:red"]),
    legend=False,
    ax=ax[1],
)

for a in ax.flat:
    a.axhline(0, color="darkslategrey", linewidth=1)
    a.tick_params(axis="x", labelsize=11)
ax[0].tick_params(axis="y", labelsize=11)
ax[0].axvline("85", color="darkslategrey", linewidth=1)
ax[1].axvline("120", color="darkslategrey", linewidth=1)
ax[0].set_xlabel("Cavern diameter [m]")
ax[1].set_xlabel("Cavern height [m]")
ax[0].set_ylabel("Difference in storage capacity [%]")
ax[0].legend(title=None)

sns.despine()
plt.tight_layout()
plt.savefig(
    os.path.join("graphics", "fig_sensitivity.jpg"),
    format="jpg",
    dpi=600,
)
plt.show()
