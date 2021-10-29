# hikari
hikari is a simple Python3.6+ package for manipulating basic crystallographic
files: mainly .hkl, but also .res, .cif and, by extension, .fcf.

The following section contains brief explanation of how to install 
and use hikari. For full description please see the documentation.

## Getting started

Currently, hikari is not registered in PyPI python package.
In order to register hikari in your python3 namespace,
firstly download or clone the repository and create new virtual environment.

    $ mkvirtualenv -p /usr/bin/python3.6 hikari-venv

While inside your virtual environment and in the repository folder run:

    $ pip install -r requirements.txt
    
to install the dependencies. After running python from the virtual environment,
the `hikari` module should be available in the namespace via `import hikari`.

## Usage

For the sake of usage, hikari can be essencially understood as divided
into two sub-modules: dataframes and scripts.
Dataframes contain object responsible for basic manipulation of files,
for example the `hikari.dataframes.CifFrame` is responsible for
reading, modifying and writing the Crystal Information Files.
Scripts, on the other hands, contain ready to use sets of dataframe
instructions and aim to solve a complete problem, like reformatting the file
or evaluating data completeness in certain experiment.

## Author

This software is made by Daniel Tchoń, dtchon@chem.uw.edu.pl, and distributed
under MIT license. Dr hab. Anna Makal and prof. dr hab. Krzysztof Woźniak are
acknowledged for helpful discussions about crystallography behind the subject.
Dr Jarosław Kalinowski and mgr inż. Damian Tchoń are acknowledged
for insightful tips on project structure and code optimisation.
