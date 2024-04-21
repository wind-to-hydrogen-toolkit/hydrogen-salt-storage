# Storage of hydrogen from offshore wind in salt caverns

Optimising production and long-term bulk storage of hydrogen from offshore wind in salt caverns in the Irish Sea.

## Acknowledgements

The [H-Wind](https://www.marei.ie/project/h-wind) project was supported by a research grant from [Science Foundation Ireland (SFI)](https://www.sfi.ie/) under Grant No. 12/RC/2302 – P2 and by the industry consortium members: [DP Energy](https://dpenergy.com/), [Equinor](https://www.equinor.com/), [ESB](https://esb.ie/), and [Gas Networks Ireland](https://www.gasnetworks.ie/).

## Installation

This project uses [Python](https://www.python.org/) ≥ 3.11.

### Installation of the project (including notebooks) from source

Clone this Git repository:

```sh
git clone https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage.git
cd hydrogen-salt-storage
```

Create a virtual environment and install all requirements:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

To run tests:

```sh
python -m pytest --cov
```

### Installation of the `h2ss` module

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage/archive/refs/heads/main.zip
```

## Documentation

Documentation is available at: <https://wind-to-hydrogen-toolkit.gitlab.io/hydrogen-salt-storage>.

To build the documentation locally:

```sh
cd docs && make html
```

To clean build the documentation locally:

```sh
cd docs && make clean html
```

## Licence

Copyright 2023-2024 N. Streethran

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

<https://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
