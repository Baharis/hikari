# hikari

[![PyPI version](https://badge.fury.io/py/hikari-toolkit.svg)](https://badge.fury.io/py/hikari-toolkit)
[![codecov](https://codecov.io/gh/Baharis/hikari/branch/master/graph/badge.svg?token=SWKKW0LSKQ)](https://codecov.io/gh/Baharis/hikari)
[![CodeFactor](https://www.codefactor.io/repository/github/baharis/hikari/badge)](https://www.codefactor.io/repository/github/baharis/hikari)
[![Documentation Status](https://readthedocs.org/projects/hikari/badge/?version=stable)](https://hikari.readthedocs.io/en/stable/?badge=stable)

hikari is a simple Python3.6+ package for manipulating basic crystallographic
files: mainly .hkl, but also .res, .cif and, by extension, .fcf.

The following section contains brief explanation of how to install 
and use hikari. For full description please see the documentation.

## Getting started

Hikari is registered in PyPI under the name `hikari-toolkit`.
In order to start working with the package, simply install it using:

    $ pip install hikari-toolkit

Since it runs on python version 3.6 or newer and requires specific versions
of some popular packages such as numpy, you might be interested in
using hikari in a virtual environment, which might be created using:

    $ mkvirtualenv -p /usr/bin/python3.6 hikari-venv

After running python from this virtual environment,
the package should be available in the namespace via `import hikari`.

## Usage

For the sake of usage, hikari is essencially divided into a few sub-modules,
including dataframes, symmetry, utility, and scripts.
Dataframes contain object responsible for basic manipulation of files,
for example the `hikari.dataframes.CifFrame` is responsible for
reading, modifying and writing the Crystal Information Files.
Scripts, on the other hands, contain ready to use sets of dataframe
instructions and aim to solve a complete problem, like reformatting the file
or evaluating data completeness in certain experiment.
In the majority of cases, you should be more interested in the latter.

## Author

This software is made by Daniel Tchoń, dtchon@chem.uw.edu.pl, and distributed
under MIT license. Dr hab. Anna Makal and prof. dr hab. Krzysztof Woźniak are
acknowledged for helpful discussions about crystallography behind the subject.
Dr Jarosław Kalinowski and mgr inż. Damian Tchoń are acknowledged
for insightful tips on project structure and code optimisation.
