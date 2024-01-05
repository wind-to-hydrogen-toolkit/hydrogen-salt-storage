# Hydrogen storage in salt caverns in the Irish Sea

## Notebooks

Notebooks can be viewed [here](https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/tree/ipynb/notebooks).

1. Kish Basin halite: [data][nbkbdata], [stats][nbkbstats]
1. [Proposed offshore wind farms][nbowf]
1. [Land boundary and basemaps][nbboundary]

### Constraints data

1. [Exploration wells][nbwells]
1. [Frequently used shipping routes][nbshipping]
1. [Shipwrecks][nbshipwrecks]
1. [Subsea cables][nbcables]

### Cavern capacity calculation

1. [Cavern generation in zone of interest][nbcaverns]
1. [Caverns with constraints applied][nbconstrained]
1. [Cavern storage capacity][nbcapacity]

### Optimisation and analysis

1. [Weibull parameters of wind speeds][nbweibull]
1. [Optimisation][nboptimisation]
1. [Sensitivity analysis][nbsensitivity]

## Datasets

Dataset | Source | Metadata | Format | Version
-- | -- | -- | -- | --
Wind Farms (Foreshore Process) | Department of Housing, Local Government and Heritage | [Link][owf]; [Link][owf2] | Shapefile | 2023-09-11
Exploration Wells in the Irish Offshore | Petroleum Affairs Division / Department of the Environment, Climate and Communications | [Link][wells]; [Link][wells2] | Shapefile | 2021-02-01
Shipwrecks in Irish Waters | INFOMAR / Geological Survey Ireland / Marine Institute / National Monuments Service / Department of Housing, Local Government and Heritage | [Link][shipwrecks] | Shapefile | 2023-10-24
Frequently Used Routes (300 gross tonnes and above) | Department of Housing, Local Government and Heritage / EMODnet | [Link][shippingroutes]; [Link][shippingroutes2] | Shapefile | 2022-07-19
Irish Sea Subsea Cables and Renewable Energy Structures | KIS-ORCA | [Link][cables] | Olex Fishing Plotter File | 2023-02-27
Weibull Parameters of Wind Speeds 2001 to 2010 - 150m above ground level | Sustainable Energy Authority of Ireland | [Link][weibull]; [Link][weibull2] | Shapefile | 2021-06-30
Provinces - National Statutory Boundaries - 2019 - Ungeneralised | Ordnance Survey Ireland / Tailte Éireann | [Link][boundary]; [Link][boundary2] | Shapefile | 2022-07-07
Kish Basin Halite | HYSS | [Link][hyss] | DAT | 2023-07
Bathymetry | EMODnet | [Link][emodnet] | GeoTIFF | 2022

## Installation

This project uses Python 3.11.

Create a virtual environment and install all requirements:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

To run tests:

```sh
python -m pytest
```

To generate a coverage report with the tests:

```sh
coverage run -m pytest && coverage report -m
```

[owf]: https://data.gov.ie/dataset/wind-farms-foreshore-process
[owf2]: https://data-housinggovie.opendata.arcgis.com/maps/housinggovie::wind-farms-foreshore-process
[wells]: https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2171
[wells2]: https://data.gov.ie/dataset/exploration-wells-in-the-irish-offshore
[hyss]: https://hyss.ie
[emodnet]: https://emodnet.ec.europa.eu/en/bathymetry
[shipwrecks]: https://isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.5131
[shippingroutes]: https://data.gov.ie/dataset/frequently-used-routes-300-gross-tonnes-and-above1
[shippingroutes2]: https://data-housinggovie.opendata.arcgis.com/maps/housinggovie::frequently-used-routes-300-gross-tonnes-and-above
[cables]: https://kis-orca.org/downloads/
[boundary]: https://data.gov.ie/dataset/provinces-osi-national-statutory-boundaries-2019
[boundary2]: https://data-osi.opendata.arcgis.com/maps/osi::provinces-national-statutory-boundaries-2019
[weibull]: https://data.gov.ie/dataset/weibull-parameters-wind-speeds-2001-to-2010-150m-above-ground-level
[weibull2]: https://gis.seai.ie/wind/
[nbkbdata]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_dat_files.ipynb
[nbkbstats]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_stats.ipynb
[nbowf]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/offshore_wind_farms.ipynb
[nbwells]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/offshore_exploration_wells.ipynb
[nbweibull]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/weibull_parameters_wind_speeds.ipynb
[nbshipping]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/frequent_shipping_routes.ipynb
[nbshipwrecks]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/shipwrecks.ipynb
[nbcables]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/subsea_cables.ipynb
[nboptimisation]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/wind_farm_optimisation.ipynb
[nbsensitivity]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_sensitivity.ipynb
[nbcaverns]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_caverns.ipynb
[nbconstrained]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_caverns_constrained.ipynb
[nbcapacity]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_caverns_capacity.ipynb
[nbboundary]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/ireland_boundary_basemaps.ipynb
