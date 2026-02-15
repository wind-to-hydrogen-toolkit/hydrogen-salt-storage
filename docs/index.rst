.. hydrogen-salt-storage documentation master file, created by
   sphinx-quickstart on Sun Dec 31 21:35:00 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Optimising production and long-term bulk storage of hydrogen from offshore wind in salt caverns in the Irish Sea
================================================================================================================

✨ Streethran, N., Byrne, K., White, J., O'Neill, N., and Leahy, P. (2026). 'Optimising production and long-term bulk storage of hydrogen from offshore wind in subsurface salt caverns', *Journal of Energy Storage*, 153, p. 120988. Available at: https://doi.org/10.1016/j.est.2026.120988.

Source code is available on GitHub: https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage.

Acknowledgements
----------------

This research was carried out as part of the `H-Wind <https://www.marei.ie/project/h-wind/>`_ and `HYSS <https://hyss.ie/>`_ projects.
The H-Wind project was supported by a research grant from `Taighde Éireann – Research Ireland <https://www.researchireland.ie/>`_ under Grant No. 12/RC/2302 – P2 and by the industry consortium members: `DP Energy <https://dpenergy.com/>`_, `Equinor <https://www.equinor.com/>`_, `ESB <https://esb.ie/>`_, and `Gas Networks Ireland <https://www.gasnetworks.ie/>`_.
The HYSS project was supported by a grant from `Sustainable Energy Authority of Ireland (SEAI) <https://www.seai.ie/>`_ and `Geological Survey Ireland <https://www.gsi.ie/>`_ under the SEAI Research, Development & Demonstration Funding Programme 2021, `Grant No. 21/RDD/725 <https://www.seai.ie/seai-research/research-database/research-projects/details/hydrogen-salt-storage-assessment-hyss>`_.

Contents
--------

.. toctree::
   :maxdepth: 2

   datasets
   methods
   notebooks
   map

Installation
------------

This project uses `Python <https://www.python.org/>`_ ≥ 3.11.

Clone the Git repository:

.. code-block:: shell

   git clone https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage.git
   cd hydrogen-salt-storage

Create a virtual environment and install all requirements:

.. code-block:: shell

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip setuptools wheel
   python -m pip install -r requirements.txt

Run tests:

.. code-block:: shell

   python -m pytest --cov

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
