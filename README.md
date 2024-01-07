# Storage of hydrogen from offshore wind in salt caverns

Optimising production and long-term bulk storage of hydrogen from offshore wind in salt caverns in the Irish Sea

This research was supported by a research grant from Science Foundation Ireland (SFI) under Grant No. 12/RC/2302 – P2 and by the H-Wind academic-industry consortium members: DP Energy, ESB, Equinor, and Gas Networks Ireland.

## Documentation

Documentation is available at: <https://hydrogen_salt_storage.readthedocs.io>

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

## Licence

Copyright 2023-2024 N. Streethran

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

<https://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
