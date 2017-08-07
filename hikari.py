"""Set of .hkl reflection file utility tools

Usage:
  hikari draw2d <hkl_path> [-i=format] [-r=path | -c=path] [-s=name]
            [--itosigma] [--plane=plane] [--reduce=times] [--scale=scale]
  hikari modify <hkl_path> [-i=format] [-o=format] [--sep] [-r=path | -c=path]
            [--dac=angle] [--pos=pos] [--thin=falloff] [--trim=res]
  hikari -h | --help

Options:
  -c=path           Path of .CIF file with cell and lambda data [default: None].
  -h --help         Show this help screen.
  -i=format         Format of input .HKL file     [default: (4, 4, 4, 8, 8, 4)].
  -o=format         Format of output .HKL file    [default: (4, 4, 4, 8, 8, 4)].
  -r=path           Path of .INS/.RES file with cell and lambda [default: None].
  -s=name           Name of saved .png and .hkl, without ext. [default: hikari].
  --dac=angle       Opening angle of the DAC in degrees, if any [default: None].
  --itosigma        Visualise I/s(I) of reflections using transparency.
  --plane=plane     Miller indices of drawn plane      [default: ('h', 'k', 0)].
  --pos=pos         Crystal plane adjacent to DAC diamond  [default: (1, 0, 0)].
  --reduce=reduce   Multiplicity of same-reflections averaging     [default: 0].
  --thin=falloff    Random reflection exclusion falloff factor  [default: None].
  --trim=res        High resolution reflection cutoff, in A^-1  [default: None].
  --scale=scale     Scale factor for reflections' size           [default: 1.0].
  --sep             Separate individual columns of hkl using a whitespace.
"""

import docopt
import numpy as np
from bin import hikari
from bin import readers

options = docopt.docopt(__doc__)

# SETTING READ AND WRITE FORMAT
for index, key in enumerate(hikari.Pattern.reading_parameters):
    hikari.Pattern.reading_parameters[key] = eval(options['-i'])[index]
for index, key in enumerate(hikari.Pattern.writing_parameters):
    hikari.Pattern.writing_parameters[key] = eval(options['-o'])[index]
hikari.Pattern.writing_separator = options['--sep']

# GENERATING NECESSARY OBJECTS
c = hikari.Crystal()
p = hikari.Pattern()
p.read(options['<hkl_path>'])

# READING CRYSTAL INFORMATION FROM RES
if options['-r'] != 'None':
    res = readers.read_res(options['-r'])
    p.edit_wavelength(res['CELL'][0])
    c.a_d = res['CELL'][1]
    c.b_d = res['CELL'][2]
    c.c_d = res['CELL'][3]
    c.al_d = np.radians(res['CELL'][4])
    c.be_d = np.radians(res['CELL'][5])
    c.ga_d = np.radians(res['CELL'][6])

# READING CRYSTAL INFORMATION FROM CIF
elif options['-c'] != 'None':
    cif = readers.read_res(options['-r'])
    p.edit_wavelength(str(cif['_diffrn_radiation_type']).strip('\',"'))
    c.a_d = readers.ustrip(cif['_cell_length_a'])
    c.b_d = readers.ustrip(cif['_cell_length_b'])
    c.c_d = readers.ustrip(cif['_cell_length_c'])
    c.al_d = np.radians(readers.ustrip(cif['_cell_angle_alpha']))
    c.be_d = np.radians(readers.ustrip(cif['_cell_angle_beta']))
    c.ga_d = np.radians(readers.ustrip(cif['_cell_angle_gamma']))

if options['draw2d']:
    p.place(c)
    for times in range(eval(options['--reduce'])):
        p.reduce()
    p.draw(color_scheme='gist_rainbow', dpi=600, itosigma=options['--itosigma'],
         legend=True, projection=eval(options['--plane']), savefig=True,
         savename=options['-s'], scale=eval(options['--scale']), showfig=None)

elif options['modify']:
    p.place(c)
    if options['--trim'] != 'None':
        p.trim(eval(options['--trim']))
    if options['--thin'] != 'None':
        p.thin_out(eval(options['--thin']))
    if options['--dac'] != 'None':
        p.dac(eval(options['--dac']), eval(options['--pos']))
    p.write(options['-s']+'.hkl')
