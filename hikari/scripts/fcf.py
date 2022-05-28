import numpy as np
from matplotlib import pyplot, cm
from scipy.optimize import minimize
from scipy.special import erfinv
from scipy.stats import norm

from hikari.dataframes import HklFrame, ResFrame
from hikari.symmetry import SG
from hikari.utility import make_abspath


def baycon_plot(x_key='ze', y_key='si',
                a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0,
                input_path='shelx.fcf',
                input_format='shelx_fcf',
                input_wavelength='MoKa',
                output_path='baycon.png'):
    """
    For a given .fcf file prepare a bayesian conditional probability plot
    between x_key and y_key.

    :param x_key: Parameter of HklFrame which will be placed on x axis
    :type x_key: str
    :param y_key: Parameter of HklFrame which will be placed on x axis
    :type y_key: str
    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param input_path: Path to the input .fcf file.
    :type input_path: str
    :param input_format: Format of the input .fcf file. For reference see
        :meth:`hikari.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :param output_path: Path to the output .png file.
    :type output_path: str
    """
    no_of_bins = 10
    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(make_abspath(input_path), input_format)
    p.place()
    p.calculate_fcf_statistics()
    x = p.table.loc[:, x_key].rank(pct=True).to_numpy()
    y = p.table.loc[:, y_key].rank(pct=True).to_numpy()
    bins = np.zeros(shape=(no_of_bins, no_of_bins))
    lims = [-1.e-8] + [(i + 1) / no_of_bins for i in range(no_of_bins)]
    for i in range(no_of_bins):
        for j in range(no_of_bins):
            bins[i, j] = ((lims[i] < x) & (x <= lims[i+1]) &
                          (lims[j] < y) & (y <= lims[j+1])).sum()
    n_avg = len(x) / no_of_bins ** 2
    chi2 = np.sum((bins - n_avg) ** 2 / n_avg)
    fig = pyplot.figure()
    ax = fig.add_subplot(111, aspect='equal')
    pyplot.xlim(0, 1)
    pyplot.ylim(0, 1)
    h = ax.hist2d(x, y, bins=no_of_bins, alpha=0.25, cmap=cm.get_cmap('PiYG'))
    cb = pyplot.colorbar(h[3], ax=ax)
    cb.set_label('Number of observations')
    ax.scatter(x=x, y=y, s=5.0, c='#000080', marker='.', alpha=0.75)
    pyplot.title('Bayesian CoNditional probability, chi2 = {:.2f}'.format(chi2))
    pyplot.xlabel('"' + x_key + '" rank')
    pyplot.ylabel('"' + y_key + '" rank')
    pyplot.tight_layout()
    pyplot.savefig(fname=make_abspath(output_path), dpi=300)


def observed_vs_calculated_plot(input_path='shelx.fcf',
                                input_format='shelx_fcf',
                                output_path='Io_vs_Ic.png'):
    p = HklFrame()
    p.read(make_abspath(input_path), input_format)
    icalc = p.table.loc[:, 'Ic'].to_numpy()
    iobs = p.table.loc[:, 'I'].to_numpy()
    i_min = min(np.min(icalc[icalc > 0]), np.min(iobs[iobs > 0]))
    i_max = max(np.max(icalc[icalc > 0]), np.max(iobs[iobs > 0]))
    fig = pyplot.figure()
    ax = fig.add_subplot(111)  # , aspect='equal'
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([i_min, i_max])
    ax.set_ylim([i_min, i_max])
    ax.plot(np.linspace(0, i_max), np.linspace(0, i_max), '-k', lw=1, zorder=0)
    ax.scatter(x=icalc, y=iobs, s=5.0, c='r', marker='.', alpha=0.75, zorder=10)
    pyplot.title('Calculated vs observed intensities plot')
    pyplot.xlabel('I_cal')
    pyplot.ylabel('I_obs')
    pyplot.tight_layout()
    pyplot.savefig(fname=make_abspath(output_path), dpi=300)


def normal_probability_plot(input_path='shelx.fcf',
                            input_format='shelx_fcf',
                            output_path='Io_vs_Ic.png'):

    # scale factors
    a = 0.1000
    b = 0.0

    p = HklFrame()
    p.read(make_abspath(input_path), input_format)
    i_obs = p.table.loc[:, 'I'].to_numpy()
    i_calc = p.table.loc[:, 'Ic'].to_numpy()
    si = p.table.loc[:, 'si'].to_numpy()
    p = 1/3 * i_obs + 2/3 * i_calc
    si = np.sqrt(si ** 2 + (a * p) ** 2 + b * p)

    # expected delta m
    def delta_m(f1, f2, k, si1, si2):
        return np.sort((f1 - k * f2) / np.sqrt(si1 ** 2 + k **2 * si2 ** 2))

    def sum_of_delta_m_squared(k):
        return np.sum(delta_m(i_obs, i_calc, k, si, np.zeros_like(si)) ** 2)

    def scale_factor():
        return minimize(sum_of_delta_m_squared, x0=np.array([1.0])).x[0]

    experiment_delta_m = delta_m(f1=i_obs, f2=i_calc, k=scale_factor(),
                               si1=si, si2=np.zeros_like(si))
    experiment_delta_m = experiment_delta_m / np.std(experiment_delta_m)

    # simulated delta m
    uniform = (np.arange(len(experiment_delta_m))+0.5) / len(experiment_delta_m)
    simulated_delta_m = [erfinv(-1 + 2 * q) for q in uniform]

    # drawing the plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    pyplot.hist(experiment_delta_m, bins=100, density=True)
    ax.scatter(experiment_delta_m, simulated_delta_m, s=5.0, c='r', marker='.',
               alpha=0.75, zorder=10)
    ax.plot(np.linspace(-3, 3), np.linspace(-3, 3), '-k', lw=1, zorder=0)
    pyplot.plot(6 * uniform - 3, norm.pdf(6 * uniform - 3))
    pyplot.title('npp')
    pyplot.xlabel('delta_m experiment')
    pyplot.ylabel('delta_m simulated')
    pyplot.tight_layout()
    pyplot.savefig(fname=make_abspath(output_path), dpi=300)


def fcf_descriptors(input_path='shelx.fcf', input_format='shelx_fcf'):
    # scale factors
    a = 0.1000
    b = 0.0

    p = HklFrame()
    p.read(make_abspath(input_path), input_format)
    i_obs = p.table.loc[:, 'I'].to_numpy()
    i_calc = p.table.loc[:, 'Ic'].to_numpy()
    si = p.table.loc[:, 'si'].to_numpy()
    p = 1/3 * i_obs + 2/3 * i_calc
    si_weighted = np.sqrt(si ** 2 + (a * p) ** 2 + b * p)
    ze = (i_obs - i_calc) / si_weighted
    f_calc = np.sqrt(np.abs(i_calc)) * np.sign(i_calc)
    f_obs = np.sqrt(np.abs(i_obs)) * np.sign(i_obs)
    one_over_sf = (2 * abs(i_obs) ** 0.5) / si

    r1 = np.sum(np.abs(f_obs - f_calc)) / np.sum(np.abs(f_obs))
    wr2 = np.sqrt(
        np.sum(np.abs(si_weighted * np.abs(i_obs - i_calc) ** 2)) /
        np.sum(np.abs(si_weighted * i_obs ** 2)))
    awr2 = np.sqrt(
        (np.mean((i_obs - i_calc) ** 2) / np.mean(si_weighted ** 2)) /
        np.mean((i_obs / si_weighted) ** 2))
    gof_if_alpha_equal_one = np.sqrt(np.mean(ze ** 2))
    agof_if_alpha_equal_one = np.sqrt(
        np.mean((i_obs - i_calc) ** 2) /
        np.mean(si_weighted ** 2))

    print('R1    = {:f}'.format(r1))
    print('wR2   = {:f}'.format(wr2))
    print('awR2  = {:f}'.format(awr2))
    print('GoF*  = {:f}'.format(gof_if_alpha_equal_one))
    print('aGoF* = {:f}'.format(agof_if_alpha_equal_one))


def calculate_sample_form_factors(a, b, c, al, be, ga, space_group, res_path):
    """
    Estimate and print selected IAM XRD form factors for given crystal structure

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param space_group: Short Hermann-Mauguin name or index of space group.
        For details see table in hikari.symmetry.space_groups.
    :type space_group: str or int
    :param res_path: Absolute or relative path to the input .res file.
    :type res_path: str
    :return: None
    :rtype: None
    """
    r = ResFrame()
    r.read(make_abspath(res_path))
    r.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    hkl = np.array([(0, 0, 0), (1, 1, 1), (2, 2, 2), (2, 0, 0), (0, 0, 3),
                (1, 0, 1), (1, 1, 8), (5, 0, 2), (4, 4, 0), (2, 0, 6),
                (2, 0, 1), (2, 0, 2), (2, 0, 3), (2, 0, 4), (2, 0, 5),
                (5, 9, 9), (0, 0, 10), (0, 2, 10), (0, 4, 10)])
    f = r.form_factor(np.array(hkl), SG[space_group])
    f2 = f * np.conj(f)
    for _hkl, _f, _f2 in zip(hkl, f, f2):
        print(f'{_hkl}: {_f2:12f} --- {_f}')


if __name__ == '__main__':
    # calculate_sample_form_factors(a=5.64109, b=5.64109, c=5.64109,
    #                               al=90, be=90, ga=90, space_group='Fm-3m',
    #                               res_path='~/x/NaCl/cifmaking/NaCl_more_res.res')
    calculate_sample_form_factors(a=7.210241, b=16.487567, c=11.279203,
                                  al=90, be=90, ga=90, space_group='Pnma',
                                  res_path='~/x/HP/2oAP/_/_.res')
