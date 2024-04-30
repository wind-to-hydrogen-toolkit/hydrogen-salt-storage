#!/usr/bin/env python
# coding: utf-8

# # Net-to-gross
#
# Data from HYSS for Kish Basin (<https://hyss.ie/>)

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

data = pd.DataFrame(
    {
        "gross": [190, 201, 219, 144, 278, 271, 120],
        "NTG": [0.68, 0.34, 0.39, 0.41, 0.41, 0.51, 0.36],
        "well": [
            "33/21-1",
            "33/21-1",
            "33/21-1",
            "33/21-1",
            "33/17-2A",
            "33/17-2A",
            "33/17-1",
        ],
        "halite": [
            "Preesall",
            "Mythop",
            "Rossall",
            "Fylde",
            "Rossall",
            "Fylde",
            "Fylde",
        ],
    }
)

data.sort_values(by=["gross", "NTG"], inplace=True)

data

model = sm.OLS(data["NTG"], sm.add_constant(data["gross"]))
results = model.fit()

print(results.summary())

b, m = results.params
r = results.rsquared

g = sns.lmplot(data=data, x="gross", y="NTG")
plt.text(215, 0.65, f"$y = {m:.5f}x {b:+.5f}$\n$R^2 = {r:.5f}$", fontsize=11.5)
g.set_axis_labels("Gross halite thickness [m]", "Net-to-gross ratio")
plt.show()

g = sns.lmplot(data=data, x="gross", y="NTG", hue="halite")
# g.set(xlim=(0, 500), ylim=(0, 1))
g.set_axis_labels("Gross halite thickness [m]", "Net-to-gross ratio")
plt.show()

data.describe()

# ## Fylde only

fylde = data[data["halite"] == "Fylde"]

model = sm.OLS(fylde["NTG"], sm.add_constant(fylde["gross"]))
results = model.fit()

print(results.summary())

b, m = results.params
r = results.rsquared

g = sns.lmplot(data=fylde, x="gross", y="NTG")
plt.text(210, 0.75, f"$y = {m:.5f}x {b:+.5f}$\n$R^2 = {r:.5f}$", fontsize=11.5)
g.set_axis_labels("Gross halite thickness [m]", "Net-to-gross ratio")
plt.show()

data.describe()

# ## Linear regression


def net_to_gross(gross):
    y = m * gross + b
    return min(y, 0.75)


ntg = []
gross_thickness = np.arange(0, 700, step=1)

for x in gross_thickness:
    ntg.append(net_to_gross(x))

df = pd.DataFrame({"gross": gross_thickness, "NTG": ntg})

df["net"] = df["gross"] * df["NTG"]

df.describe()

net_to_gross(1000)

print(f"{net_to_gross(300):.5f}")

ax = sns.scatterplot(
    data=data,
    x="gross",
    y="NTG",
    hue="halite",
    zorder=3,
    palette="mako",
    s=75,
)
df.plot(
    x="gross",
    y="NTG",
    zorder=1,
    color="slategrey",
    linestyle="dashed",
    label=f"$y = \\min({m:.5f}x {b:+.5f}, 0.75)$\n$R^2 = {r:.5f}$",
    linewidth=2,
    ax=ax,
)
ax.set_xlabel("Gross halite thickness [m]")
ax.set_ylabel("Net-to-gross ratio")
sns.despine()
ax.set(xlim=(0, 700), ylim=(0, 0.8))
ax.xaxis.grid(True, linewidth=0.25)
ax.yaxis.grid(True, linewidth=0.25)
plt.legend(title=None, loc="lower right")
plt.tight_layout()
# plt.savefig(
#     os.path.join("graphics", "fig_net_to_gross.jpg"),
#     format="jpg",
#     dpi=600,
# )
plt.show()
