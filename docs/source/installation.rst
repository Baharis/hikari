Installation
====================

Introduction
********************

Kesshou is a simple Python3 package which aims to
read, edit, analyse and visualise various crystallographic files.
Since it's beginning it has been developed mainly as a self-teaching project,
but due to its simplicity, adaptability, and relative accessibility
it can offer a lot of functionality for common users.
At the moment the library is in an early stage of development, so placement,
naming and functionality of various modules or functions is a subject to change.

At it's current stage of development the software is focused mainly on
modifying and visualizing the .hkl single crystal reflection files.
It's basic capabilities include:

- Reading and writing various formats of .hkl files,
- Calculating simple .hkl statistics,
- Visualising the .hkl contents,
- Removing data from existing .hkl files,
- Predicting quality of high-pressure experiments

In order to fully utilize library's functionality, you will need:
- Python 3.4 or newer,
- Gnuplot
- CCDC Mercury
- Elementary knowledge of Python.

While both *Gnuplot* and *Mercury* are not strictly necessary,
they are used to provide visualisation of quality higher than Python tools.
Some of the kesshou methods output `hkl.res` and `.gnu` files,
which can be then imported to *Mercury* or *Gnuplot*
in order to improve the data visualisation quality.

Download
********************

Kesshou is yet to be registered as a public python package, so its installation
requires a few additional steps to register it in Python's namespace.
Please mind that the software has been designed on Linux Mate for Python 3.5
and thus the instruction will reference to this exact operating system
and this exact version of Python. Both Windows and Mac operating systems
have analogous tools to perform necessary, presented steps.

Start by cloning or downloading Kesshou to the directory of your choice.
By default, the new *project* directory will be called "Kesshou",
and it will contain a folder with an actual *package* called "kesshou".

Virtual environment
********************

It is recommended to use Kesshou exclusively inside a virtual environment
such as `virtualenvwrapper <http://virtualenvwrapper.readthedocs.io>`_
in order to avoid version conflicts of with system installed packages.
This will be a subject to change when the software will be morphed into
a python package, but at the moment it protects the system from undesired
conflicts. New virtual environment on Linux can be created using:

.. code-block:: bash

    $ mkvirtualenv -p /usr/bin/python3.5 kesshou_venv

Where `/usr/bin/python3.5` is a path to the python3 executable.

Virtualenvwrapper is available on Mac OS. It can be also found on Windows
under the name virtualenvwrapper-win. Once the virtual environment is
created, you will be automatically placed inside, which in majority of terminals
should be signified using additional braces at the beginning of the prompt:

.. code-block:: bash

    (kesshou-venv) $

You can now freely leave and re-enter the newly-created environment
using `workon kesshou_venv` and `deactivate`, respectively.

Dependencies
********************

The software depends on a few popular external python modules, such as
`numpy <https://numpy.org>`_, `pandas <https://pandas.pydata.org/>`_ or
`matplotlib <https://matplotlib.org/>`_,
all of which are very consistent, developed and popular.
A full list of direct dependencies can be checked by investigating the
`requirements.txt` file in the project directory.

In order to install the dependencies, first and foremost please make sure
that you are working inside the virtual environment. Please navigate
to the "Kesshou" folder and install all dependencies using:

.. code-block:: bash

    (kesshou-venv) $ pip install -r requirements.txt

Pip is a package manager for Python modules (packages)
and is available on all three mentioned operating systems.
The installation process might take a few minutes depending on the
availability of the libraries on your local machine and internet connection.

Finally it is necessary to tell your Python executable where to look for
the library itself. While in your virtual environment,
append the "Kesshou" *project* directory to the PYTHONPATH variable using:

.. code-block:: bash

    (kesshou-venv) $ export PYTHONPATH='/absolute/path/to/folder/Kesshou/'

At this point kesshou should be available just like any Python package
whenever you run Python from the virtual environment. While the package does
not yet have test module, you can check if the istallation have been successful
by running Python and executing:

.. code-block:: python

    import kesshou

While you might happen to see some *hopefully irrelevant* warnings,
a lack of :code:`ImportError` should signify that the installation
have been performed correctly.
