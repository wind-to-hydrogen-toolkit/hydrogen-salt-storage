#!/usr/bin/env python
# coding: utf-8

# # Sensitivity analysis

# In[ ]:


import os
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from h2ss import compare

# In[ ]:


cavern_diameter = np.arange(45, 106, step=5)
cavern_height = np.arange(80, 321, step=20)


# In[ ]:


def generate_sensitivity_data(diameter, height):
    """Generate data to perform sensitivity analysis"""
    os.makedirs(os.path.join("data", "sensitivity"), exist_ok=True)
    ds, extent, exclusions = compare.load_all_data()
    for d, h in product(diameter, height):
        filename = os.path.join(
            "data", "sensitivity", f"sensitivity_d{d}_h{h}.csv"
        )
        if not os.path.isfile(filename):
            df, _ = compare.capacity_function(ds, extent, exclusions, d, h)
            df = df[["cavern_diameter", "cavern_height", "capacity"]]
            df.to_csv(filename)
            print(f"{filename} done!")
        else:
            print(f"{filename} exists!")


# In[ ]:


generate_sensitivity_data(cavern_diameter, cavern_height)


# In[ ]:


filelist = []
for cd, ch in product(cavern_diameter, cavern_height):
    filelist.append(
        os.path.join("data", "sensitivity", f"sensitivity_d{cd}_h{ch}.csv")
    )


# In[ ]:


cdf = pd.concat((pd.read_csv(f) for f in filelist), ignore_index=True)


# In[ ]:


cdf.drop(columns=["Unnamed: 0"], inplace=True)


# In[ ]:


cdf["cavern_height"] = cdf["cavern_height"].astype(int)


# In[ ]:


cdf.describe()


# In[ ]:


# cavern height-to-diameter-ratio
pd.Series((cdf["cavern_height"] / cdf["cavern_diameter"]).unique()).describe()


# In[ ]:


ax = sns.scatterplot(
    data=cdf,
    hue="cavern_diameter",
    y="capacity",
    x="cavern_height",
    linewidth=0,
    alpha=0.75,
    s=5,
)
sns.despine()
plt.show()


# In[ ]:


len(cdf["cavern_diameter"].unique()) == len(cdf["cavern_height"].unique())


# ## Number of caverns and total capacity

# In[ ]:


f, ax = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

data = (
    cdf.groupby(["cavern_height", "cavern_diameter"])
    .count()
    .reset_index()
    .pivot(index="cavern_height", columns="cavern_diameter", values="capacity")
    .sort_index(ascending=False)
)
sns.heatmap(
    data,
    ax=ax[0],
    cmap="mako_r",
    cbar=False,
    square=True,
    annot=True,
    fmt=",d",
    robust=True,
)

data = cdf.copy()
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

# plt.savefig(
#     os.path.join("graphics", "fig_sensitivity_heatmap.jpg"),
#     format="jpg",
#     dpi=600,
# )
plt.show()


# ## Mean capacity

# In[ ]:


data = (
    cdf.groupby(["cavern_height", "cavern_diameter"])
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

# In[ ]:


base = cdf[
    (cdf["cavern_diameter"] == 85) & (cdf["cavern_height"] == 120)
].reset_index(drop=True)


# In[ ]:


base.describe()[["capacity"]]


# In[ ]:


base_mean = base[["capacity"]].mean().values[0]


# In[ ]:


base_sum = base[["capacity"]].sum().values[0]


# In[ ]:


base_count = base[["capacity"]].count().values[0]


# In[ ]:


print(f"{base_mean:.3f}")


# In[ ]:


print(f"{base_sum:.3f}")


# In[ ]:


base_count


# ## Base diameter, varying height

# In[ ]:


dd = cdf[(cdf["cavern_diameter"] == 85)].reset_index(drop=True)


# In[ ]:


dd_mean = (
    pd.DataFrame(dd.groupby("cavern_height").mean()["capacity"] - base_mean)
    / base_mean
    * 100
).reset_index()


# In[ ]:


dd_sum = (
    pd.DataFrame(dd.groupby("cavern_height").sum()["capacity"] - base_sum)
    / base_sum
    * 100
).reset_index()


# In[ ]:


dd_count = (
    pd.DataFrame(dd.groupby("cavern_height").count()["capacity"] - base_count)
    / base_count
    * 100
).reset_index()


# ## Base height, varying diameter

# In[ ]:


dh = cdf[(cdf["cavern_height"] == 120)].reset_index(drop=True)


# In[ ]:


dh_mean = (
    pd.DataFrame(dh.groupby("cavern_diameter").mean()["capacity"] - base_mean)
    / base_mean
    * 100
).reset_index()


# In[ ]:


dh_sum = (
    pd.DataFrame(dh.groupby("cavern_diameter").sum()["capacity"] - base_sum)
    / base_sum
    * 100
).reset_index()


# In[ ]:


dh_count = (
    pd.DataFrame(
        dh.groupby("cavern_diameter").count()["capacity"] - base_count
    )
    / base_count
    * 100
).reset_index()


# ## Combined plots

# In[ ]:


dh_sum["type"] = "Total capacity"
dh_count["type"] = "Number of caverns"
dd_sum["type"] = "Total capacity"
dd_count["type"] = "Number of caverns"


# In[ ]:


f, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
sns.barplot(
    data=pd.concat([dh_count, dh_sum]),
    hue="type",
    y="capacity",
    x="cavern_diameter",
    palette=sns.color_palette(["tab:red", "tab:blue"]),
    legend=False,
    ax=ax[0],
)
sns.barplot(
    data=pd.concat([dd_count, dd_sum]),
    hue="type",
    y="capacity",
    x="cavern_height",
    palette=sns.color_palette(["tab:red", "tab:blue"]),
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
ax[1].legend(title=None, fontsize=12)
ax[1].tick_params(axis="y", left=False)

sns.despine()
plt.tight_layout()
# plt.savefig(
#     os.path.join("graphics", "Figure_9.png"),
#     format="png",
#     dpi=600,
# )
plt.show()


# In[ ]:


def colour_label(dataset):
    "Use different colours for negative and positive values"
    conditions = [(dataset["capacity"] < 0), (dataset["capacity"] >= 0)]
    choices = ["N", "P"]
    dataset["colour"] = np.select(conditions, choices)
    return dataset


# In[ ]:


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
