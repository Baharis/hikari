# Kesshou
Simple Python3 package for manipulating crystallographic files.

## Introduction

Kesshou is a simple Python3 package which aims to
read, edit, analyse and visualise various crystallographic files.
Since it's beginning it has been developed mainly as a self-teaching project,
but due to its simplicity, adaptability, and relative accessibility
it came to offer a lot of functionality for common users.
At the moment the library is in an early stage of development,
so placement and functionality of various modules or functions
is still a subject to change.
While the major bugs have been eradicated,
an unfortunate combination of variables can lead to 

At it's current stage of development
the software is focused mainly on
modifying and visualizing the .hkl single crystal reflection files.
It's basic capabilities include:
- Reading and writing various formats of .hkl files,
- Calculating simple .hkl statistics,
- Visualising the .hkl contents,
- Removing data from existing .hkl files,
- Predicting quality of high-pressure experiments

## Getting Started

### Prerequisites

In order to fully utilize library's functionality, you will need the following:
- Python 3.4 or newer,
- Gnuplot
- CCDC Mercury 
- Elementary knowledge of Python.
 
While both *Gnuplot* and *Mercury* are not strictly necessary,
they are used to provide visualisation of quality higher than Python tools.
Some of the kesshou methods output `hkl.res` and `.gnu` files,
which can be then imported to *Mercury* or *Gnuplot*
in order to improve the data visualisation quality.

### Instalation

Kesshou is yet to be registered as a public python package,
so its installation requires a few additional steps.
Please mind that the software has been designed on Linux Mate for Python 3.5
and thus the instruction will reference to this exact operating system
and this exact version of Python.

Start by cloning or downloading Kesshou to the directory of your choice.
By default, the new directory will be called "Kesshou",
and it will contain a folder with an actual package called "kesshou".
It is recommended to use Kesshou inside a virtual environment
such as [virtualenvwrapper](http://virtualenvwrapper.readthedocs.io)
in order to avoid version conflicts of with system installed packages.
Virtual environment on Linux can be created using:

    $ mkvirtualenv -p /usr/bin/python3.5 kesshou
    
In order to re-enter or leave the newly-created environment
use `workon kesshou` and `deactivate`, respectively.
**While in the kesshou virtual environment**,
please navigate to the "Kesshou" folder and install all dependencies using

    $ pip install -r requirements.txt

This installation process might take a few minutes.
Finally it is necessary to tell your Python executable where to look for
the library itself. Next, while in your virtual environment,
append the "Kesshou" directory to the PYTHONPATH variable using:

    $ export PYTHONPATH='/absolute/path/to/folder/Kesshou/'

The aforementioned steps can be performed on
both **Windows** and **Mac** using analogous tools.

### Testing the functionality

At this point you should be able too use kesshou just like any Python package.
In order to check the installation, 
run your Python executable and try importing the package using:
    
```python
import kesshou
```

You might see some irrelevant warnings which should not impend programs
functionality. If you will not see the 
`ImportError: No module named kesshou` line,
your installation has most likely been performed correctly.
In the future it is planned to include the tests to check it.


## Mini-tutorial

This mini-tutorial will revolve about reflection data gathered for
doxycycline hydrochloride and published
[here](https://doi.org/10.1515/zkri-2018-2058).
Basic crystallographic data necessary for this task
has been presented in a table below.
If you prefer to work with your own data, feel free to do so.

| Docycycline hydrochloride||
|---------------|------------|
| Space group   | P2<sub>1</sub>2<sub>1</sub>2<sub>1</sub>  |
| a /Angstrom   | 11.079     |
| b /Angstrom   | 12.674     |
| c /Angstrom   | 16.677     |
| alpha /degree | 90.000     |
| beta  /degree | 90.000     |
| gamma /degree | 90.000     |
| wavelength    | 0.71073    |

### Read and write

After activating the virtual environment, run python in command prompt
or interactive shell.
Since we are interested in the .hkl reflection file,
we will load a dataframe made specifically to work with these files.
Then we will create an object called `my_hkl`,
which will be a direct representation of our hkl file:

```pythonstub
from kesshou.dataframes.hkl import HklFrame
my_hkl = HklFrame()
```
This `HklFrame` is one of many data frames and
it contains all methods used to manipulate a single-crystal reflection files.
If, at any point, you are interested with existing frames or their methods,
do not hesitate to open a relevant file and check out its documentation.

For now our hkl data frame is completely empty - 
and we would like to fill it with some sample data.
Firstly we should *always* provide the basic information
about unit cell and wavelength, which is necessary to position reflections
in the reciprocal space. You can do this using:

```pythonstub
my_hkl.crystal.edit_cell(a=11.079, b=12.674, c=16.677, al=90, be=90, ga=90)
my_hkl.edit_wavelength(0.71073)
```

All kesshou functions accept by default Angstrom and degrees as their units
of measurement, and use two-letter abbreviations for greek letters.
Since we have provided basic information about the unit cell,
we are ready to read the .hkl file.
Since it is written in hkl4 format, we need to specify:

```pythonstub
my_hkl.read(hkl_path='example/doxy.hkl', hkl_format=4)
```

Note that your `hkl_path` might be different,
depending on where did you run your python executable.
Now the frame is filled with data, which you can check using `len(my_hkl)`.
We would now like to write this data in different format,
for example to use it in *XD* software. In this case just use:

```pythonstub
my_hkl.write(hkl_path='example/xd.hkl', hkl_format='XD')
```

You can exploit this functionality to re-format any number of hkl files
to any defined (and even custom) format.
Remember to consult either the documentation, the code or its author
when in doubt.

### Analyse and visualise

Reading and writing can be already pretty useful,
but let's say we would like to know what exactly is in our file.
Firstly, we can try to get some basic statistics using
```pythonstub
my_hkl.make_stats()
```
This function will present us with some basic statistics about
redundancy, completeness etc. in different resolution shells.
Note, that the resolution is given as distance from (000),
so multiplying it by 0.5 leads to crystallographic resolution in sin(th)/la.
As you can see, the completeness in all shells is equal roughly 25%.
This is caused by the fact that we are yet to supply the symmetry information.
In order to do that, we will firstly import the definition of all point groups,
and then pass the relevant group to the `make_stats` function:

```pythonstub
from kesshou.symmetry.pointgroup import *
my_hkl.make_stats(point_group=PG222)
```
The nomenclature behind point groups is fairly straightforward
and has been presented in *docs/pointgroup_names.txt*. New values of
completeness are much closer to 100%. Missing data comes from the fact,
that we have not considered systematic extinctions,
a subject which will not be tackled in this tutorial.

In order to actually *see* your data, it is firstly advised to 
import the style file *example/hkl.msd* into your *Mercury* software,
turn it on, and then restart the program to fully update the style.
Kesshou produces large *.res* files and then passes reflections as atoms,
so the style file makes it easier for both program and user
to handle the task at hand.
Nonetheless the computer might stutter while loading this many "atoms",
so we will also limit the data to the resolution of 0.837A:

```pythonstub
my_hkl.trim(limit=1 / 0.837)
my_hkl.to_hklres(path='example/hkl.res')
```

This will produce a new "hkl.res" file which, after some loading time,
should be shown in *Mercury*. Basic settings of visualise will
use atom size to present the I/sigma value. If the size of points is not to
your liking, you might try increasing the intensity a few times (eg. 5) using
`my_hkl.rescale('I', 5)` and running the `to_hklres` method again.

### Show dac completeness map

Basic functions present in `kesshou.dataframes` can be very powerful,
but the user should not be forced to learn programming
in order to evaluate his data. Accessing more powerful tools prepared beforehand
can be performed using "taskmasters" - the scripts, which utilize basic 
methods of dataframes to directly provide user desired data.
In this example, we will try to look for optimal Diamond Anvil Cell (dac)
orientation to pack out crystal, so we will:

```pythonstub
from kesshou.taskmasters.hkl import completeness_map
```
By going to the relevant file `kesshou/taskmasters/hkl/completeness_map`
you can see that the script is fairly long,
but in principle it does only use a bunch of methods from the dataframe.
Since each of these methods requires different variables to be provided,
this taskmaster requires a lot of input:

```pythonstub
from kesshou.taskmasters.hkl import completeness_map
completeness_map(a=11.079,
                 b=12.674,
                 c=16.677,
                 al=90,
                 be=90,
                 ga=90,
                 laue_group=PGmmm,
                 extinctions=('h00: h=2n', '0k0: k=2n', '00l: l=2n'),
                 opening_angle=35,
                 output_directory='example/',
                 output_name='doxy_cplt_map',
                 output_quality=3,
                 resolution=0.83,
                 wavelength='MoKa')
```

This code should run a few minutes,
so we have some time to talk about its functionality.
It will assign the unit cell of *a, b, c, al, be, ga*
and a *wavelength* to our data. For a series of different directions
(whose number depends on *output_quality*) it will fill the hkl dataframe
with all reflections which should theoretically fit inside a DAC
rotated in this direction (excluding extincted ones).
The (single, not double) *opening_angle* and *resolution* of data
can be additionally controlled.
Then the script will check for an actual completeness in provided "laue_group"
and will output the result to the .lst and .dat files.
Using these file, a 2D map in .png format will be created to give user
an insight about how the crystal should be placed.
A *.gnu* file will also be proveded,
which should produce a map of slightly better quality.

## Author

This software is made by Daniel Tchoń, dtchon@chem.uw.edu.pl.
The author thanks dr hab. Anna Makal and prof. dr hab. Krzysztof Woźniak
for helpful discussion about crystallography behind the subject.
The author also thanks dr Jarosław Kalinowski and mgr inż. Damian Tchoń
for helpful discussion about Python code optimisation.
