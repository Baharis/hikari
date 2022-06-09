# hikari

[![PyPI version](https://img.shields.io/pypi/v/hikari-toolkit)](https://pypi.org/project/hikari-toolkit/)
[![Python version](https://img.shields.io/pypi/pyversions/hikari-toolkit.svg)](https://www.python.org/downloads/release/python-3615/)
[![codecov](https://codecov.io/gh/Baharis/hikari/branch/master/graph/badge.svg?token=SWKKW0LSKQ)](https://codecov.io/gh/Baharis/hikari)
[![CodeFactor](https://www.codefactor.io/repository/github/baharis/hikari/badge)](https://www.codefactor.io/repository/github/baharis/hikari)
[![Documentation Status](https://readthedocs.org/projects/hikari/badge/?version=stable)](https://hikari.readthedocs.io/en/stable/?badge=stable)
[![tests](https://github.com/Baharis/hikari/actions/workflows/codecov.yml/badge.svg?branch=master)](https://github.com/Baharis/hikari/actions/workflows/codecov.yml)

hikari is a simple Python3.6+ package intended for manipulating and running 
scripts on basic crystallographic files:
.hkl, .fcf, .cif, and to some extent .res and .lst.

The following section contains a brief explanation of how to install 
and use hikari. For a full description please see
[the documentation](https://hikari.readthedocs.io/en/stable/?badge=stable).

## Getting started

Hikari is registered in PyPI under the name `hikari-toolkit`.
In order to start working with the package, install it using:

    $ pip install hikari-toolkit

Since it runs on Python 3.6+ and requires specific versions of some popular
packages such as `numpy`, you might be interested in using hikari
in a virtual environment. On Linux, it can be created using `virtualenvwrapper`:

    $ mkvirtualenv -p /usr/bin/python3.6 hikari-venv

After running python from this virtual environment,
the package should be available in the namespace via `import hikari`.

## Usage

For the sake of usage, hikari is essentially divided into a few sub-modules,
including dataframes, symmetry, utility, and scripts.
Dataframes contain objects responsible for basic manipulation of files,
for example, the `hikari.dataframes.CifFrame` is responsible for
reading, modifying and writing the Crystal Information Files.
Scripts, on the other hand, contain ready to use sets of dataframe
instructions and aim to solve certain problems, like reformatting the file
or evaluating data completeness in an experiment. In the majority of cases,
you will be most likely more interested in the latter.

## Author

This software is made by
[Daniel Tchoń](https://www.researchgate.net/profile/Daniel-Tchon),
and distributed under an MIT license. It is in constant development and all
tips, suggestions, or contributions are welcome and can be sent
[here](mailto:dtchon@chem.uw.edu.pl).
If you have utilised `hikari` in academic work, please consider citing 
[this article](https://doi.org/10.1107/S2052252521009532).

Dr hab. Anna Makal and prof. dr hab. Krzysztof Woźniak are
acknowledged for helpful discussions about crystallography behind the subject.
Dr Jarosław Kalinowski and mgr inż. Damian Tchoń are acknowledged
for insightful tips on project structure and code optimisation.
