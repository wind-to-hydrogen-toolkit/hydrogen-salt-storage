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

### Cavern generation and capacity calculation

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
Wind Farms (Foreshore Process) | Department of Housing, Local Government and Heritage | [Link][owf] | Shapefile | 2023-09-11
Exploration Wells in the Irish Offshore | Petroleum Affairs Division / Department of the Environment, Climate and Communications | [Link][wells]; [Link][wells2] | Shapefile | 2021-02-01
Shipwrecks in Irish Waters | INFOMAR / National Monuments Service / Department of Housing, Local Government and Heritage | [Link][shipwrecks] | Shapefile | 2023-10-24
Frequently Used Routes (300 gross tonnes and above) | Department of Housing, Local Government and Heritage / EMODnet | [Link][shippingroutes] | Shapefile | 2022-07-19
Weibull Parameters of Wind Speeds 2001 to 2010 - 150m above ground level | Sustainable Energy Authority of Ireland | [Link][weibull] | Shapefile | 2021-06-30
Provinces - OSi National Statutory Boundaries - 2019 - Ungeneralised | Ordnance Survey Ireland / Tailte Éireann | [Link][boundary] | Shapefile | 2022-07-07
Kish Basin Halite | HYSS | [Link][hyss] | DAT | 2023-07

[owf]: https://data.gov.ie/dataset/wind-farms-foreshore-process
[wells]: https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2171
[wells2]: https://data.gov.ie/dataset/exploration-wells-in-the-irish-offshore
[hyss]: https://hyss.ie
[shipwrecks]: https://isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.5131
[shippingroutes]: https://data.gov.ie/dataset/frequently-used-routes-300-gross-tonnes-and-above1
[boundary]: https://data.gov.ie/dataset/provinces-osi-national-statutory-boundaries-2019
[weibull]: https://data.gov.ie/dataset/weibull-parameters-wind-speeds-2001-to-2010-150m-above-ground-level
[nbkbdata]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_dat_files.ipynb
[nbkbstats]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_stats.ipynb
[nbowf]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/wind_farms_foreshore_process.ipynb
[nbwells]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/offshore_exploration_wells.ipynb
[nbweibull]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/weibull_parameters_wind_speeds.ipynb
[nbshipping]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/frequent_shipping_routes.ipynb
[nbshipwrecks]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/heritage_sites.ipynb
[nbcables]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/subsea_cables.ipynb
[nboptimisation]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/wind_farm_optimisation.ipynb
[nbsensitivity]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_sensitivity.ipynb
[nbcaverns]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_caverns.ipynb
[nbconstrained]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_caverns_constrained.ipynb
[nbcapacity]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/kish_basin_caverns_capacity.ipynb
[nbboundary]: https://github.com/nmstreethran/hydrogen-salt-storage-optimisation/blob/ipynb/notebooks/ireland_boundary_basemaps.ipynb
