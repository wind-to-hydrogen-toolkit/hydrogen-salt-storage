#!/usr/bin/env python
# coding: utf-8

# # Sensitivity analysis

import os
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from h2ss import compare

cavern_diameter = np.arange(45, 106, step=5)
cavern_height = np.arange(80, 321, step=20)


def generate_sensitivity_data(cavern_diameter, cavern_height):
    """Generate data to perform sensitivity analysis"""
    os.makedirs(os.path.join("data", "sensitivity"), exist_ok=True)
    ds, extent, exclusions = compare.load_all_data()
    for d, h in product(cavern_diameter, cavern_height):
        filename = os.path.join(
            "data", "sensitivity", f"sensitivity_d{d}_h{h}.csv"
        )
        if not os.path.isfile(filename):
            df = compare.capacity_function(ds, extent, exclusions, d, h)
            df.to_csv(filename)
            print(f"{filename} done!")
        else:
            print(f"{filename} exists!")


# generate_sensitivity_data(cavern_diameter, cavern_height)

filelist = []
for d, h in product(cavern_diameter, cavern_height):
    filelist.append(
        os.path.join("data", "sensitivity", f"sensitivity_d{d}_h{h}.csv")
    )

df = pd.concat((pd.read_csv(f) for f in filelist), ignore_index=True)

df.drop(columns=["Unnamed: 0"], inplace=True)

df["cavern_height"] = df["cavern_height"].astype(int)

df.describe()

# cavern height-to-diameter-ratio
pd.Series((df["cavern_height"] / df["cavern_diameter"]).unique()).describe()

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

len(df["cavern_diameter"].unique()) == len(df["cavern_height"].unique())

# ## Number of caverns and total capacity

f, ax = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

data = (
    df.groupby(["cavern_height", "cavern_diameter"])
    .count()
    .reset_index()
    .pivot(index="cavern_height", columns="cavern_diameter", values="capacity")
    .sort_index(ascending=False)
)
sns.heatmap(
    data,
    ax=ax[0],
    cmap="rocket_r",
    cbar=False,
    square=True,
    annot=True,
    fmt=",d",
    robust=True,
)

data = df.copy()
data["capacity"] = data["capacity"] / 1000
data = (
    data.groupby(["cavern_height", "cavern_diameter"])
    .sum()
    .reset_index()
    .pivot(index="cavern_height", columns="cavern_diameter", values="capacity")
    .sort_index(ascending=False)
)
sns.heatmap(
    data,
    ax=ax[1],
    cmap="mako_r",
    cbar=False,
    square=True,
    annot=True,
    robust=True,
    fmt=".2f",
)

for a in ax.flat:
    a.tick_params(axis="y", labelsize=11.5, left=False)
    a.tick_params(axis="x", labelsize=11.5, bottom=False)
    a.set_xlabel("Cavern diameter [m]")

ax[0].set_ylabel("Cavern height [m]")
ax[1].set_ylabel("")
ax[0].set_title("Number of caverns")
ax[1].set_title("Total capacity [TWh]")
plt.tight_layout()

plt.savefig(
    os.path.join("graphics", "fig_sensitivity_heatmap.jpg"),
    format="jpg",
    dpi=600,
)
plt.show()

# ## Mean capacity

data = (
    df.groupby(["cavern_height", "cavern_diameter"])
    .mean()
    .reset_index()
    .pivot(index="cavern_height", columns="cavern_diameter", values="capacity")
    .sort_index(ascending=False)
)
f, ax = plt.subplots(figsize=(7, 7))
sns.heatmap(
    data,
    ax=ax,
    cmap="mako_r",
    # cbar_kws={"label": "Mean capacity [GWh]"},
    cbar=False,
    square=True,
    annot=True,
    robust=True,
    fmt=".2f",
)
ax.set_xlabel("Cavern diameter [m]")
ax.set_ylabel("Cavern height [m]")
ax.tick_params(axis="y", labelsize=11.5, left=False)
ax.tick_params(axis="x", labelsize=11.5, bottom=False)
plt.title("Mean capacity [GWh]")
# ax.plot(8.5, 10.5, "*", markersize=10, color="black")
plt.tight_layout()
plt.show()

# ## Base case

base = df[
    (df["cavern_diameter"] == 85) & (df["cavern_height"] == 120)
].reset_index(drop=True)

base.describe()[["capacity"]]

base_mean = base[["capacity"]].mean().values[0]

base_sum = base[["capacity"]].sum().values[0]

base_count = base[["capacity"]].count().values[0]

print(f"{base_mean:.3f}")

print(f"{base_sum:.3f}")

base_count

# ## Base diameter, varying height

dd = df[(df["cavern_diameter"] == 85)].reset_index(drop=True)

dd_mean = (
    pd.DataFrame(dd.groupby("cavern_height").mean()["capacity"] - base_mean)
    / base_mean
    * 100
).reset_index()

dd_sum = (
    pd.DataFrame(dd.groupby("cavern_height").sum()["capacity"] - base_sum)
    / base_sum
    * 100
).reset_index()

dd_count = (
    pd.DataFrame(dd.groupby("cavern_height").count()["capacity"] - base_count)
    / base_count
    * 100
).reset_index()

# ## Base height, varying diameter

dh = df[(df["cavern_height"] == 120)].reset_index(drop=True)

dh_mean = (
    pd.DataFrame(dh.groupby("cavern_diameter").mean()["capacity"] - base_mean)
    / base_mean
    * 100
).reset_index()

dh_sum = (
    pd.DataFrame(dh.groupby("cavern_diameter").sum()["capacity"] - base_sum)
    / base_sum
    * 100
).reset_index()

dh_count = (
    pd.DataFrame(
        dh.groupby("cavern_diameter").count()["capacity"] - base_count
    )
    / base_count
    * 100
).reset_index()

# ## Combined plots

dh_sum["type"] = "Total capacity"
dh_count["type"] = "Number of caverns"
dd_sum["type"] = "Total capacity"
dd_count["type"] = "Number of caverns"

f, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
sns.barplot(
    data=pd.concat([dh_count, dh_sum]),
    hue="type",
    y="capacity",
    x="cavern_diameter",
    palette=sns.color_palette(["tab:red", "tab:blue"]),
    ax=ax[0],
)
sns.barplot(
    data=pd.concat([dd_count, dd_sum]),
    hue="type",
    y="capacity",
    x="cavern_height",
    palette=sns.color_palette(["tab:red", "tab:blue"]),
    legend=False,
    ax=ax[1],
)
for a in ax.flat:
    a.axhline(0, color="darkslategrey", linewidth=1)
    a.yaxis.grid(True, linewidth=0.25)
    # a.set(ylim=(-100, 200))
ax[0].axvline("85", color="darkslategrey", linewidth=1, linestyle="dashed")
ax[1].axvline("120", color="darkslategrey", linewidth=1, linestyle="dashed")
ax[0].set_xlabel("Cavern diameter [m]")
ax[1].set_xlabel("Cavern height [m]")
ax[0].set_ylabel("Difference [%]")
ax[0].legend(title=None, fontsize=11)
ax[1].tick_params(axis="y", left=False)

sns.despine()
plt.tight_layout()
plt.savefig(
    os.path.join("graphics", "fig_sensitivity.jpg"),
    format="jpg",
    dpi=600,
)
plt.show()


def colour_label(df):
    conditions = [(df["capacity"] < 0), (df["capacity"] >= 0)]
    choices = ["N", "P"]
    df["colour"] = np.select(conditions, choices)
    return df


f, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
sns.barplot(
    data=colour_label(dh_mean),
    y="capacity",
    x="cavern_diameter",
    hue="colour",
    palette=sns.color_palette(["tab:red", "tab:blue"]),
    legend=False,
    ax=ax[0],
)
sns.barplot(
    data=colour_label(dd_mean),
    y="capacity",
    x="cavern_height",
    hue="colour",
    palette=sns.color_palette(["tab:red", "tab:blue"]),
    legend=False,
    ax=ax[1],
)

for a in ax.flat:
    a.axhline(0, color="darkslategrey", linewidth=1)
    a.yaxis.grid(True, linewidth=0.25)
ax[0].axvline("85", color="darkslategrey", linewidth=1, linestyle="dashed")
ax[1].axvline("120", color="darkslategrey", linewidth=1, linestyle="dashed")
ax[0].set_xlabel("Cavern diameter [m]")
ax[1].set_xlabel("Cavern height [m]")
ax[0].set_ylabel("Difference in mean storage capacity [%]")
ax[1].tick_params(axis="y", left=False)

sns.despine()
plt.tight_layout()
plt.show()
