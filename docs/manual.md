<pre>
  _  __             _                 
 | |/ /___  ___ ___| |__   ___  _   _ 
 | ' // _ \/ __/ __| '_ \ / _ \| | | |
 | . \  __/\__ \__ \ | | | (_) | |_| |
 |_|\_\___||___/___/_| |_|\___/ \__,_|
</pre>

## Introduction
Kesshou is a simple Python3 package which aims to
read, edit, analyse and visualise various crystallographic files.
Since it's beginning it have been developed mainly as a self-teaching project,
but due to its simplicity, adaptability and relative accessibility
it came to offer a lot of options for common users.

At it's current stage of development the software is focused mainly on
modifying the .hkl files. It's basic capabilities include:
- Reading and writing various formats of .hkl:
    - popular formats such as hkl4 or hkl5
    - program-specific formats such as xd or tonto
    - custom formats using user-made definition
- Calculating simple .hkl statistics:
    - I/sigma, redundancy, completeness vs resolution
    - assuming both no symmetry and given point group
    - completeness in DAC conditions
- Visualising the .hkl contents:
    - in 2D, using slicing and plotting scatter plots 
    - in 3D, using mercury as external visualiser
- Removing data from existing .hkl files:
    - By 'trimming' the file to given resolution
    - By randomly deleting certain percent of data
    - By simulating Diamond Anvil Cell shadow

## Object-oriented strategy
The package has been written using object-oriented programming 
design philosophy. While as a self-teaching project it is certainly flawed,
it bases on some simple and easy-to-understand principles. 
Firstly, each file type is handled using a related **dataframe**.
In principle all information regarding .hkl files in stored and manipulated
using a hkl dataframe, or more precisely a `HklFrame` stored in `dataframe.hkl`
directory.

The `HklFrame` contains all methods used to manipulate a reflection files.
In general, the methods are in-place, which means they do not return anything
and modify themselves instead.
In order to obtain a new `HklFrame` object
one has to create a clear copy using `HklFrame()`
or copy existing one using `deepcopy` from `copy` package.
The following list contains some of the more useful and handy methods
written for `HklFrame`:
- `HklFrame.read()` reads and interprets file `hkl_path`, 
whose format is `hkl_format`;
- `HklFrame.write()` outputs 
the file to `hkl_path` in `hkl_format`, using blank `columns_separator` if needed;
- `HklFrame.reduce()` performs a very simple procedure of
merging reflections with common h, k, l integers;
- `HklFrame.place()` positions all signals in reciprocal space
which is necessary for some other methods;
- `HklFrame.dac()` removes data which should be shadowed
by a DAC of given `opening_angle` and orientation;
- `HklFrame.make_stats()` returns basic statistics in resolution `bins` 
using given `point_group`.

Since learning all the methods and writing the code from scratch
might be tiresome for some users, 
package `Kesshou` is bundled with a folder `scripts`
which contains some basic routines ready to analyze the data.
Those scripts should suffice for majority of potential software applications,
and some of them will be explained in the following chapters.

## Reading and writing

The most basic operation which any user would perform on the data is 
importing and exporting the reflection information in and out of the program.
The script containing only those basic procedures,
called `hkl_change-format.py`, 
will be of use to learn the very basics of Kesshou.
The script is divided into three separate segments,
which will be discussed below.

Firstly, the file contains the import statements.
In this example only one import statement is necessary,
as all the methods we are planning to use belong to the `HklFrame` class.
Therefore, this part consists of ony one uncommented line:
`from kesshou.dataframes.hkl import HklFrame`.
Since in all discussed examples we will work with hkl files,
this statement will be necessary in almost every analysed script.

The second part of the file consists of the 'variables' field.
In order to increase the overall clarity of the code,
this area has been separated from the main code,
and should be essentially the only place for user input.
The type of input which is required from the user
should be clearly explained in the comments of each file.
In this case, the script asks for four different variables:
- `input_hkl_path` - a path to the .hkl file which is to be loaded;
- `output_hkl_path` - a path where the exported .hkl file will be saved;
- `input_hkl_format` - format of the loaded .hkl file
- `output_hkl_format` - format of the exported .hkl file

Finally, the third part of the file consists of raw code, which executes
thanks to the 'import statements' and according to the 'variables'.
In the case of `hkl_change-format.py` essentially
only three operations are performed:
1) `p = HklFrame()` - an HklFrame object to store the data called 'p' is created;
1) `p.read(input_hkl_path, input_hkl_format)` - data from `input_hkl_path`
 is imported into the 'p';
1) `p.write(output_hkl_path, output_hkl_format)` - data from 'p'
is exported into `output_hkl_path`.

In order to execute this (and every other) script, 
it is suggested to run the script either directly from a python editor
or in a command line (using `python3 hkl_change-format.py`),
since at this point majority of scripts do not produce readable output
and dump the results to the command line instead. 

The main functionality of this very script is to re-format any hkl file.
For example, tonto software cannot work with standard hkl4 and instead
requires user to prepare a very specific type of reflection file.
Automating this process saves users their valuable time and can be
quickly executed in Kesshou.

As a final note in this section please note that
the scripts have been prepared as a kind of a guideline and 
the software does not perform any quality checks on input variables e.g,:
it will not check whether input text or number makes sense.
In order to make your input is reasonable,
please see where in the code is it used and refer to the documentation
of each individual function itself.

## Trimming data to DAC limits

Now let's move to some more sophisticated instructions.
The second script we are going to tackle,
`dac_cut-hkl-to-DAC-limit-and-resolution-limit.py`,
does much more than just reformatting a file;
Nonetheless its structure is identical to the first one's.
It starts from the import statements, which this time apart from
the `HklFrame` import also `numpy` and `copy` - 
two common libraries used by Kesshou to create and copy objects.

```
from kesshou.dataframes.hkl import HklFrame
import numpy as np
from copy import deepcopy
```

The script then asks user to input a lot of information concerning
crystal structure and desired output. These information,
and the reasoning behind a need for them, include:
- unit cell parameters - 
    to correctly position reflections in reciprocal space
- orientation matrix or DAC-perpendicular vector -
    to correctly position DAC relative to crystal
- values of opening angles - 
    for each of them a separate set of files will be prepared
- resolution limit to be imposed on each dataset as desired
- input and output file details

Please note that this time in input you are additionally asked to provide
the wavelength, as it is necessary to determine a shape of DAC shadow.
Moreover, instead of output file name you are asked to provide a
directory and base name, as the script will prepare numerous files
based on this input and provide name for each of them would be more tiring.

The variables field in this case might look something like this:

```
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 5.0857
unit_cell_b = 11.804
unit_cell_c = 5.4606
unit_cell_al = 90
unit_cell_be = 111.98
unit_cell_ga = 90

# Crystal orientation matrix from .cif file
UB_11 = -0.0483900000
UB_12 = 0.0237766000
UB_13 = -0.0253660000
UB_21 = 0.0588330000
UB_22 = 0.0310240000
UB_23 = -0.0129997000
UB_31 = -0.0119742000
UB_32 = -0.0528324000
UB_33 = -0.0431912000

# OR DAC perpendicular vector, if better suited
v1 = 0.0
v2 = 1.0
v3 = 0.0
use_vector_instead_of_orientation_matrix = True

# Opening angle in degrees
pressure_cell_oa = [35, ]

# Resolution limit as sin(th)/la, A-1
reslim = 0.8

# Input details
input_hkl_path = '/input/directory/glycine.hkl'
input_hkl_format = 'xd'
input_hkl_wavelength = 'MoKa'

# Output details
output_name = 'glycine'
output_directory = '/output/directory/glycine.hkl'
output_hkl_format = 4
```

This input will cause the script to import 
the hkl data from `input_hkl_path` file,
correctly orient the reflections in in the reiciprocal space,
cut out those which should be shadowed by Diamond Anvil Cell
(whose orientation is perpendicular to 0k0, according to [v1, v2, v3]),
remove reflections further that 1.6A-1 in Ewald construction,
and output new hkl and some simple 2D graphics of the new dataset.

```
# Prepare HklFrame object
p = HklFrame()
p.read(input_hkl_path, input_hkl_format)
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(input_hkl_wavelength)
if not use_vector_instead_of_orientation_matrix:
    p.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
                                        (UB_21, UB_22, UB_23),
                                        (UB_31, UB_32, UB_33)))
p.drop_zero()
p.place()

# Prepare list of interesting projections
projections = (('h', 'k', 0), ('h', 0, 'l'), (0, 'k', 'l'),
               ('h', 'k', 1), ('h', 1, 'l'), (1, 'k', 'l'))

# Trim data to resolution if necesary
p.trim(reslim)


# Draw projections before dac operation
q = deepcopy(p)
q.reduce()
for projection in projections:
    output_png_path = output_directory + output_name + '_full_' + \
                      str(projection[0]) + \
                      str(projection[1]) + \
                      str(projection[2])
    q.draw(colored='m', projection=projection, savepath=output_png_path)

# Cut and draw projections for subsequent DAC opening angles
try:
    _ = iter(pressure_cell_oa)
except TypeError:
    pressure_cell_oa = list(pressure_cell_oa)
for oa in pressure_cell_oa:
    q = deepcopy(p)
    if use_vector_instead_of_orientation_matrix:
        vector = np.array((v1, v2, v3))
        q.dac(opening_angle=oa, vector=vector)
    else:
        q.dac(opening_angle=oa)
    output_hkl_path = output_directory + output_name +\
                      '_oa' + str(oa)[0:6] + '.hkl'
    q.write(hkl_path=output_hkl_path, hkl_format=output_hkl_format)
    q.reduce()
    q.drop_zero()
    for projection in projections:
        output_png_path = output_directory + output_name + \
                          '_oa' + str(oa)[0:6] + '_' + \
                          str(projection[0]) + \
                          str(projection[1]) + \
                          str(projection[2])
        q.draw(colored='m', projection=projection, savepath=output_png_path)
```

## Calculating completeness in DAC