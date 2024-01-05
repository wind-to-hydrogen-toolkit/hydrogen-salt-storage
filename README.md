# Hydrogen storage in salt caverns in the Irish Sea

## Documentation

Documentation is available at: <https://hydrogen-salt-storage.readthedocs.io>

## Notebooks

Notebooks can be viewed [here](https://github.com/nmstreethran/hydrogen-salt-storage/tree/ipynb/notebooks).

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
