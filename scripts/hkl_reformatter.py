"""Changes .hkl file to a desired format

Usage:
  hkl_reformatter.py [--method=METHOD] [--plot=PLOT] [--min_snr=MIN_SNR] [--min_n=MIN_N] [--type=TYPE] [--temperature-C=TEMP_C] [--temperature-unc=TEMP_UNC] SPECTRA_FILE

Options:
  --method=METHOD             peak finding method: peakutils or find_peaks_cwt [default: find_peaks_cwt]
  --plot=PLOT                 file name for plot
  --min_snr=MIN_SNR           min_snr parameter for scipy.signal.find_peaks_cwt
  --min_n=MIN_N               minimum number of points for fitting profile
  --type=TYPE                 file type: txt, SpectraSuite, spc [default: txt]
  --temperature-C=TEMP_C      temperature in Celsius [default: 24.85]
  --temperature-unc=TEMP_UNC  temperature uncertainty [default: 0.1]
"""

import json
from os import path

from docopt import docopt

from dacxrd.algo.find_peaks import find_peaks, Method
from dacxrd.algo.j_appl_phys_110_043513_2011 import pressure as pressure_j_appl_phys_110_043513_2011
from dacxrd.readers.spectral_data import SpectralDataFileFormat, read_spectral_data
from dacxrd.spectrum import Unit


arguments = docopt(__doc__)

method = Method(arguments['--method'])

method_opts = {}
interpolate_opts = {}

if arguments['--min_snr']:
    method_opts['min_snr'] = int(arguments['--min_snr'])

if arguments['--min_n']:
    interpolate_opts['min_n'] = int(arguments['--min_n'])

T_Celsius = float(arguments['--temperature-C'])
T_unc = float(arguments['--temperature-unc'])

ftype = SpectralDataFileFormat(arguments['--type'])
fname = arguments['SPECTRA_FILE']
with open(fname, "rb") as byte_stream:
    spectral_data = read_spectral_data(byte_stream, ftype=ftype, unit=Unit.NM)

peak_R2, peak_R1 = find_peaks(
    spectral_data.data,
    method=method,
    method_opts=method_opts,
    interpolate_opts=interpolate_opts
)

peak_R2_data = peak_R2.format()
peak_R1_data = peak_R1.format()

T_Kelvin = T_Celsius + 273.15
P_j_appl_phys_110_043513_2011 = pressure_j_appl_phys_110_043513_2011(
    R1=peak_R1_data["center"],
    R1_unc=peak_R1_data["center_unc"],
    T=T_Kelvin,
    T_unc=T_unc
)

root, ext = path.splitext(path.basename(fname))
result = dict(
    file_root=root,
    file_ext=ext,
    file_dir=path.dirname(fname),
    R2=peak_R2_data,
    R1=peak_R1_data,
    T_Kelvin=T_Kelvin,
    T_unc=T_unc,
    P_j_appl_phys_110_043513_2011=P_j_appl_phys_110_043513_2011
)
if ftype in (SpectralDataFileFormat.SpectraSuite, SpectralDataFileFormat.SPC):
    result['meta'] = spectral_data.meta

result_txt = json.dumps(result, sort_keys=True, indent=4)

print(result_txt)

if arguments["--plot"]:
    import matplotlib.pyplot as plt
    from dacxrd.plot_peaks import plot_peaks


    plt.figure(1, figsize=(18, 10), dpi=100)
    plot_peaks(spectral_data.data, peak_R2, peak_R1)
    plt.savefig(arguments["--plot"], papertype='a3')

