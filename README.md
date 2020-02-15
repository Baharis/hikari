# Kesshou
kesshou [kes-show] is a simple Python3.5+ package for 
manipulating basic crystallographic files such as .hkl, .res or .cif.

The following section contains brief explanation of how to install 
and use kesshou. For full description please see the documentation.

## Getting started

At this moment python is not registered as python package for confidentiality.
In order to register kesshou in your python namespace, firstly download or clone
the repository and create new virtual environment.

    $ mkvirtualenv -p /usr/bin/python3.5 kesshou-venv

While inside your virtual environment and in the repository folder run:

    $ pip install -r requirements.txt
    
to install the dependencies. After running python from the virtual environment,
the `kesshou` module should be available in the namespace via `import Kesshou`.

## Usage

Kesshou is essencially divided into two sub-modules: dataframes and scripts.
Dataframes contain object responsible for simple manipulation of files,
for example the `kesshou.dataframes.CifFrame` is responsible for
reading, modifying and writing the Crystal Information Files.
Scripts, on the other hands, contain ready to use sets of dataframe
instructions and aim to solve a complete problem, like reformatting the file
of calculating diffraction data completeness.

## Author

This software is made by Daniel Tchoń, dtchon@chem.uw.edu.pl.
The author thanks dr hab. Anna Makal and prof. dr hab. Krzysztof Woźniak
for helpful discussion about crystallography behind the subject.
The author also thanks dr Jarosław Kalinowski and mgr inż. Damian Tchoń
for helpful discussion about code optimisation.
