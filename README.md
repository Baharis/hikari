# Kesshou
Simple file modification and visualisation toolkit for crystallography.
Available for Python 3.5+.

* module **hikari**: modification and visualisation of .hkl files
* module **pruby**: pressure calculation based on Ruby Pressure Scale spectrum

Available under the MIT License.

### Dependencies
* [docopt](http://docopt.org/)
* [matplotlib, numpy](http://www.scipy.org)
* [PeakUtils](http://pythonhosted.org/PeakUtils)
* [pyqtgraph](http://www.pyqtgraph.org/)
* [uncertainties](http://pythonhosted.org/uncertainties/)

### Installation and usage

It is recommended to use Kesshou inside a Python virtual environment
([virtualenvwrapper](http://virtualenvwrapper.readthedocs.io))
to avoid version conflicts of its dependencies with system installed packages.
In order to create virtual environment run in the command line:

    $ mkvirtualenv -p /usr/bin/python3.5 kesshou

And then install dependencies:

    $ pip install -r requirements.txt

In order to run the scripts, download or clone the repository on your local disc 
and enter the virtual environment with:

    $ workon kesshou


### Examples

For a simple .hkl file example.hkl located in test_data:

```
   1   0   0  10000.   1000.   1
   1   0   0   5000.    500.   2
   0   1   0   2500.    250.   3
  -1   0   0   1250.    125.   4
   0  -1   0   -625.    62.5   5
   0   0   1   312.5   31.25   6
   0   0   0      1.      1.   1
```

Running a command:
 
    python3 hikari.py draw2d --scale=2 test_data/example.hkl
    
Results in a visualisation of a default (hk0) plane:

![example.png](test_data/example1.png?raw=true)

Resulting image can be
easily modified using appropriate commands, for example:

    python3 hikari.py draw2d --plane="('h', 0, 'l')" --reduce=1 --scale=2 test_data/example.hkl

will visualise reflections' redundancy on (h0l) cross-section:

![example.png](test_data/example2.png?raw=true)

Reflection file may be easily modified as well. Executing

    python3 hikari.py modify -o "(2, 2, 2, 6, 6, 2)" test_data/example.hkl

Will reformat the .hkl format to the shorted form:

```
 1 0 010000.1000.0 1
 1 0 05000.0 500.0 2
 0 1 02500.0 250.0 3
-1 0 01250.0 125.0 4
 0-1 0-625.0  62.5 5
 0 0 1 312.5 31.25 6
 0 0 0   0.0   0.0 1
```

In order to see all available options execute `hikari.py --help`.
