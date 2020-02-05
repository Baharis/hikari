# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
import matplotlib.pyplot as pp
import numpy as np
import scipy.optimize as opt

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
"""This script takes a .txt spectrum file and fits to the data in file
an arbitrary, polynomial function of order below specified threshold.
Then, it yields a maximum of given function - maximum luminescence value

input_file - full or relative path to the input file
max_order - max order of fitted polynomial, greater = slower & more precise"""

# Input details
input_file = '/home/dtchon/x/HP/RFpirazB/18_08_16_RFpirazB/f4.txt'
fit_method = 'camel' #''consecutive'  # or 'camel'
max_order = 20  # irrelevant for camel fit

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #


# create generic dict-like for keeping the data
class Bunch(dict):
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self


# import the spectrum from the file and sort it
dots = Bunch()
dots.data = np.loadtxt(input_file, dtype=(float, float))
dots.data = dots.data[dots.data[:, 0].argsort()]
dots.x = dots.data[:, 0]
dots.y = dots.data[:, 1]
dots.xspan = np.max(dots.x) - np.min(dots.x)
dots.yspan = np.max(dots.y) - np.min(dots.y)

# generate sigmas - greater change means higher sigma
dots.sigmas = np.zeros(len(dots.x))
for index in range(1, len(dots.x) - 1):
    dots.sigmas[index] = abs(dots.y[index - 1] - dots.y[index]) + \
                         abs(dots.y[index + 1] - dots.y[index])
dots.sigmas[0] = dots.sigmas[1]
dots.sigmas[len(dots.sigmas) - 1] = dots.sigmas[len(dots.sigmas) - 2]
dots.sigmas_max = max(dots.sigmas)
dots.sigmas_rescaled = dots.sigmas / dots.sigmas_max
dots.importance = np.ones(len(dots.x)) - dots.sigmas_rescaled

# draw the dots
pp.minorticks_on()
pp.grid(b=True, which='major', color='gray', alpha=0.2)
pp.grid(b=True, axis='x', which='minor', color='gray', alpha=0.1)
pp.grid(b=True, axis='y', which='major', color='gray', alpha=0.1)
pp.tick_params(axis='x', which='minor', bottom='on')
colors = [(0, 0, 0, imp) for imp in dots.importance]
pp.scatter(dots.x, dots.y, marker='.', color=colors)


def consecutive_fit(dots):

    # generate consecutive polynomials to fit into the curve
    degree_range = list(range(0, max_order+1))
    curve = dict()
    prediction = [np.average(dots.y)]
    for d in degree_range:
        # fit the curve of given degree to the dots
        curve[d] = Bunch()
        curve[d].func = lambda x, *params: np.polyval(params, x)
        curve[d].results = opt.curve_fit(curve[d].func, xdata=dots.x, ydata=dots.y,
                          p0=prediction, sigma=dots.sigmas, absolute_sigma=False)
        curve[d].coefficients = curve[d].results[0]
        curve[d].curve = lambda x: np.polyval(curve[d].coefficients, x)
        curve[d].y = [curve[d].curve(x) for x in dots.x]
        prediction = np.concatenate((np.zeros(1), curve[d].coefficients))

        # score the function
        residual_list = [abs(curve[d].curve(x)-y) for x, y in zip(dots.x, dots.y)]
        curve[d].residuals = sum(residual_list)
        curve[d].score = np.dot(np.array(residual_list), dots.importance)

    # decide which function is best and output it
    scores = [curve[d].score for d in degree_range]
    best_degree = np.argmin(scores)
    curve = curve[best_degree]
    curve.degree = best_degree

    # find the maximum of best polynomial and draw it
    curve.xpeak = opt.minimize_scalar(
        lambda x: -np.polyval(curve.coefficients, x),
        bounds=(min(dots.x), max(dots.x)), method='Bounded').x
    curve.ypeak = np.polyval(curve.coefficients, curve.xpeak)
    pp.plot(curve.xpeak, curve.ypeak,
        color='green', marker='v', markersize='8',
        label='order: ' + str(curve.degree) + ', height: ' + str(curve.ypeak))

    # draw the best fit
    pp.plot(dots.x, curve.y, marker=None, color='olive', linestyle='-')
    pp.fill_between(x=dots.x, y1=0, y2=curve.y, color='olive', alpha=0.1)

    # throw everything out
    print('order: ' + str(curve.degree) + '\nheight: ' + str(curve.ypeak))

    return curve


def camel_fit(dots):

    # prepare a prediction for a dromedaries-type curve
    gauss_x1 = min(dots.x) + dots.xspan / 8
    gauss_y1 = min(dots.y) + 3 * dots.yspan / 4
    gauss_s1 = dots.xspan / 10
    gauss_x2 = min(dots.x) + dots.xspan / 2
    gauss_y2 = min(dots.y) + 1 * dots.yspan / 8
    gauss_s2 = dots.xspan / 5
    prediction = [gauss_x1, gauss_y1, gauss_s1, gauss_x2, gauss_y2, gauss_s2]
    print(prediction)

    # fit dromedaries-type curve (gaussian-gaussian + gaussian-gaussian to data
    curve = Bunch()
    curve.g = lambda x, a, m, s: a * np.exp(-(x - m) ** 2 / (2. * s ** 2))
    curve.func = lambda x, mu1, a1, si1, mu2, a2, si2: \
        curve.g(x, a1, mu1, si1) + curve.g(x, a2, mu2, si2)
    curve.results = opt.curve_fit(curve.func, xdata=dots.x, ydata=dots.y,
                    p0=prediction, sigma=dots.sigmas, absolute_sigma=False)
    curve.coefficients = curve.results[0]
    print('babayaga')
    print(curve.coefficients)
    [mu1, a1, si1, mu2, a2, si2] = curve.coefficients
    curve.curve = lambda x: curve.func(x, mu1, a1, si1, mu2, a2, si2)
    curve.y = [curve.curve(x) for x in dots.x]
    curve.degree = 'Gauss functions'

    # find the maximum 1 and draw it
    curve.xpeak = curve.coefficients[0]
    curve.ypeak = curve.coefficients[1]
    pp.plot(curve.xpeak, curve.ypeak,
        color='green', marker='v', markersize='8',
        label='Gaussian 1, x:' + str(curve.xpeak) + ', y: ' + str(curve.ypeak))

    # find the maximum 2 and draw it
    curve.xpeak = curve.coefficients[3]
    curve.ypeak = curve.coefficients[4]
    pp.plot(curve.xpeak, curve.ypeak,
        color='orange', marker='v', markersize='8',
        label='Gaussian 2, x:' + str(curve.xpeak) + ', y: ' + str(curve.ypeak))

    # draw the best fit
    pp.plot(dots.x, curve.y, marker=None, color='olive', linestyle='-')
    pp.fill_between(x=dots.x, y1=0, y2=curve.y, color='olive', alpha=0.1)

    # throw everything out
    print('order: ' + str(curve.degree) + '\nheight: ' + str(curve.ypeak))

    return curve


if fit_method == 'consecutive':
    curve = consecutive_fit(dots)
else: # fit_method == 'dromedaries':
    curve = camel_fit(dots)



pp.legend()
pp.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# TODO add documentation and scipy?
