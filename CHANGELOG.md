# Changelog

## [2024.04.1](https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage/tree/2024.04.1)

- update task and script for converting notebooks
- update documentation link - hosted on GitLab pages
- remove sidebar from datasets page in docs
- use marine area consent data for offshore wind farms instead of foreshore data
- update datasets and metadata
- add variable cavern heights based on Hystories project scenarios
- remove redundant code
- add optimisation based on 2 GW of offshore wind dedicated to hydrogen production scenario
- update docstrings and references
- use 85 m diameter and 120 m height caverns for the base case to be consistent with HYSS
- update sensitivity analysis plots and functions
- remove cavern volume correction factors - only use shape correction
- improve print statements
- calculate distance of caverns to the nearest offshore pipeline
- compare storage capacity to the projected hydrogen demand in Ireland in 2050
- improve Weibull functions
- move additional optimisation steps (e.g. AEP, number of caverns, transmission distance, LCOT) from notebook to their own functions
- add statsmodels to pyproject.toml
- use mapclassify to generate map legend and classes (quantiles)
- improve visualisations
- remove old statistical plots of Kish basin data
- update Kish basin halite net thickness facet maps
- format notebooks and scripts
