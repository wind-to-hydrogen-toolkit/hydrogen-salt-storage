.. hydrogen-salt-storage documentation master file, created by
   sphinx-quickstart on Sun Dec 31 21:35:00 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

hydrogen-salt-storage
=====================

Optimising production and long-term bulk storage of hydrogen from offshore wind in salt caverns in the Irish Sea

Source code is available on GitHub: https://github.com/nmstreethran/hydrogen-salt-storage

This research was supported by a research grant from Science Foundation Ireland (SFI) under Grant No. 12/RC/2302 – P2 and by the H-Wind academic-industry consortium members: DP Energy, ESB, Equinor, and Gas Networks Ireland.

Contents
--------

.. toctree::
   :maxdepth: 2

   datasets
   methods
   notebooks

Installation
------------

This project uses Python 3.11.

Clone the GitHub repository:

.. code-block:: shell

   git clone https://github.com/nmstreethran/hydrogen-salt-storage.git
   cd hydrogen-salt-storage

Create a virtual environment and install all requirements:

.. code-block:: shell

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install -r requirements.txt

Run tests:

.. code-block:: shell

   python -m pytest

To generate a coverage report with the tests:

.. code-block:: shell

   coverage run -m pytest && coverage report -m

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
