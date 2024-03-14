.. hydrogen-salt-storage documentation master file, created by
   sphinx-quickstart on Sun Dec 31 21:35:00 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Storage of hydrogen from offshore wind in salt caverns
======================================================

Optimising production and long-term bulk storage of hydrogen from offshore wind in salt caverns in the Irish Sea.

Source code is available on GitHub: https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage.

Acknowledgements
----------------

The `H-Wind <https://www.marei.ie/project/h-wind/>`_ project was supported by a research grant from `Science Foundation Ireland (SFI) <https://www.sfi.ie/>`_ under Grant No. 12/RC/2302 – P2 and by the industry consortium members: `DP Energy <https://dpenergy.com/>`_, `Equinor <https://www.equinor.com/>`_, `ESB <https://esb.ie/>`_, and `Gas Networks Ireland <https://www.gasnetworks.ie/>`_.

.. image:: https://raw.githubusercontent.com/wind-to-hydrogen-toolkit/.github/main/images/logos.png

Contents
--------

.. toctree::
   :maxdepth: 2

   datasets
   methods
   notebooks

Installation
------------

This project uses `Python <https://www.python.org/>`_ 3.11.

Clone the Git repository:

.. code-block:: shell

   git clone https://github.com/wind-to-hydrogen-toolkit/hydrogen-salt-storage.git
   cd hydrogen-salt-storage

Create a virtual environment and install all requirements:

.. code-block:: shell

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt

Run tests:

.. code-block:: shell

   python -m pytest --cov

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
