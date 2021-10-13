Tutorial
============

Sample
*************

As mentioned in the installation section, the package at the moment revolves
mainly around modification of the hkl reflection file. This tutorial will
revolve about data gathered for doxycycline hydrochloride, published in 2018 in
`Zeitschrift für Kristallographie <https://doi.org/10.1515/zkri-2018-2058>`_.
Basic crystallographic data necessary for this compound has been presented
in a table below. If you prefer to work with your own diffraction data,
feel free to do so, as you should not encounter significant divergences.

+-----------------------------------------------------------+
| Doxycycline hydrochloride dihydrate                       |
+---------------+-------------------------------------------+
| Space group   | |P212121|                                 |
+---------------+-------------------------------------------+
| a /Angstrom   | 11.079                                    |
+---------------+-------------------------------------------+
| b /Angstrom   | 12.674                                    |
+---------------+-------------------------------------------+
| c /Angstrom   | 16.677                                    |
+---------------+-------------------------------------------+
| alpha /degree | 90.000                                    |
+---------------+-------------------------------------------+
| beta  /degree | 90.000                                    |
+---------------+-------------------------------------------+
| gamma /degree | 90.000                                    |
+---------------+-------------------------------------------+
| wavelength    | 0.71073                                   |
+---------------+-------------------------------------------+

.. |P212121| replace:: P 2 :sub:`1` 2 :sub:`1` 2 :sub:`1`



Read and write
**************

After activating the virtual environment, run python in command prompt
or interactive shell. Since we are interested in the .hkl reflection file,
we will load a dataframe made specifically to work with these files.
Then we will create an object called `my_hkl`,
which will be a direct representation of our hkl file:

.. code-block:: python

    from kesshou.dataframes import HklFrame
    my_hkl = HklFrame()

This `HklFrame` is one of many data frames and
it contains all methods used to manipulate a single-crystal reflection files.
If, at any point, you are interested with existing frames or their methods,
do not hesitate to open a relevant file and check out its documentation.

For now our hkl data frame is completely empty -
and we would like to fill it with some sample data.
Firstly we should *always* provide the basic information
about unit cell and wavelength, which is necessary to position reflections
in the reciprocal space. You can do this using:

.. code-block:: python

    my_hkl.edit_cell(a=11.079, b=12.674, c=16.677, al=90, be=90, ga=90)
    my_hkl.la = 0.71073


All kesshou functions accept by default Angstrom and degrees as their units
of measurement, and use two-letter abbreviations for greek letters.
Since we have provided basic information about the unit cell,
we are ready to read the .hkl file.
Since it is written in hkl4 format, we need to specify:

.. code-block:: python

    my_hkl.read(hkl_path='example/doxy.hkl', hkl_format='4')


Note that your `hkl_path` might be different,
depending on where did you run your python executable.
Now the frame is filled with data, which you can check using `len(my_hkl)`.
We would now like to write this data in different format,
for example to use it in *XD* software. In this case just use:

.. code-block:: python

    my_hkl.write(hkl_path='example/xd.hkl', hkl_format='XD')


You can exploit this functionality to re-format any number of hkl files
to any defined (and even custom) format. Remember to consult the documentation,
or its author when in doubt.



Analyse and visualise
*********************

Reading and writing can be already pretty useful,
but let's say we would like to know what exactly is in our file.
Firstly, we can try to get some basic statistics using:

.. code-block:: python

    print(my_hkl.stats())

This function will present us with some basic statistics about
redundancy, completeness etc. in different resolution shells.
Note, that the resolution is given as distance from (000),
so multiplying it by 0.5 leads to crystallographic resolution in sin(th)/la.
As you can see, the completeness in all shells is equal roughly 25%.
This is caused by the fact that we are yet to supply the symmetry information.
In order to do that, we will firstly import the definition of all point groups,
and then pass the relevant group to the `make_stats` function:

.. code-block:: python

    from kesshou.symmetry import PG
    my_hkl.stats(point_group=PG['222'])

The nomenclature behind point groups is fairly straightforward
and has been presented in *docs/pointgroup_names.txt*. New values of
completeness are much closer to 100%. Missing data comes from the fact,
that we have not considered systematic extinctions,
a subject which will not be tackled in this tutorial.

In order to actually *see* your data, it is firstly advised to
import the style file *example/hkl.msd* into your *Mercury* software,
turn it on, and then restart the program to fully update the style.
The file should work correctly on Mercury 3.10+,
if you have an older version of the software just use any other style -
you will not see custom element colours and you might expect slightly lower
performance, but otherwise all the information will be preserved.
Kesshou produces large *.res* files and then passes reflections as atoms,
so the style file makes it easier for both program and user
to handle the task at hand.
Nonetheless the computer might stutter while loading this many "atoms",
so we will also limit the data to the resolution of 0.837A:

.. code-block:: python

    my_hkl.trim(limit=1 / 0.837)
    my_hkl.to_res(path='example/hkl.res')


This will produce a new "hkl.res" file which, after some loading time,
should be shown in *Mercury*. For small files it should not take more than
a few seconds. However, if the loading time exceeds 5 minutes, the file is
most likely too large to be viewed this way.



Show dac completeness map
*************************

Basic functions present in `kesshou.dataframes` can be very powerful,
but the user should not be forced to learn programming
in order to evaluate his data. Accessing more powerful tools prepared beforehand
can be performed using the scripts, which utilize basic
methods of dataframes to directly provide user desired data.
In this example, we will try to look for optimal Diamond Anvil Cell (dac)
orientation to pack out crystal, so we will:

.. code-block:: python

    from kesshou.scripts import completeness_map

By going to the relevant file `kesshou/scripts/hkl`
you can see that the script is fairly long,
but in principle it does only use a bunch of methods from the dataframe.
Since each of these methods requires different variables to be provided,
this taskmaster requires a lot of input:

.. code-block:: python

    from kesshou.scripts import completeness_map
    completeness_map(a=11.079, b=12.674, c=16.677, al=90, be=90, ga=90,
                     laue_group=PGmmm,
                     extinctions=('h00: h=2n', '0k0: k=2n', '00l: l=2n'),
                     opening_angle=35,
                     output_directory='example/',
                     output_name='doxy_cplt_map',
                     output_quality=3,
                     resolution=0.83,
                     wavelength='MoKa')

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
A *.gnu* file will also be provided,
which should produce a map of slightly better quality.
