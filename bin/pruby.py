import os
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.optimize import curve_fit
import datetime


def calculate_pressure_shift(shift, peak_number=1, temperature_ruby=293.15, temperature_sample=293.15):
    """	Input:	ruby line shift in nm (float)
                peak_number 1 for R1, 2 for R2 (int)
                temperature of reference ruby in K (float)
                temperature of sample ruby in K (float)
        Return:	pressure in the cell in GPa (float)
        Descr.: Calculate the value described above	"""
    atmospheric_pressure = 101325.
    temperature_shift_correction = \
        wavenumber_to_wavelength(ragan_r_line_estimation(temperature_sample, peak_number)) - \
        wavenumber_to_wavelength(ragan_r_line_estimation(temperature_ruby, peak_number))
    shift -= temperature_shift_correction
    pressure_in_pascal = atmospheric_pressure + float(100000000 * (2.740 * shift * 10 + 0.0009 * shift * shift * 100))
    return pressure_in_pascal / 1000000000.


def gauss(x, *parameters):
    """	Input:	x (float)
                [A (float), mu (float), sigma (float)]
        Return:	value of gaussian with desired parameters in x (float)
        Descr.: Calculate the value described above	"""
    A, mu, sigma = parameters
    return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def plot_gaussian(dots, old_maxima, coeff_first, coeff_second, range):
    """	Input:	experimental x, y spectrum data (numpy n x 2 array)
                [rough position of R2 (float), rough position of R1 (float)]
                [A (float), mu (float), sigma(float)] for R2 gaussian
                [A (float), mu (float), sigma(float)] for R1 gaussian
                half the width of data gaussians were estimated on (float)
        Return:	None
        Descr.: Plot the dots and the gaussian estimations.	"""
    x = dots[:, 0]
    y = dots[:, 1]
    xspan = np.max(x) - np.min(x)
    yspan = np.max(y) - np.min(y)
    xmin = np.min(x) - 0.05 * xspan
    xmax = np.max(x) + 0.05 * xspan
    ymin = min(np.min(y) - 0.01 * yspan, 0)
    ymax = np.max(y) + 0.20 * yspan
    x1min, x1max = old_maxima[0] - range, old_maxima[0] + range
    x2min, x2max = old_maxima[1] - range, old_maxima[1] + range
    x1 = np.linspace(x1min, x1max, int((x1max - x1min) * 100.))
    x2 = np.linspace(x2min, x2max, int((x2max - x2min) * 100.))
    pyplot.axis([xmin, xmax, ymin, ymax])
    pyplot.grid(True, color='gray', alpha=0.1)
    # pyplot.plot((coeff_first[1], coeff_first[1]), (ymin, coeff_first[0]), 'y-')
    # pyplot.plot((coeff_second[1], coeff_second[1]), (ymin, coeff_second[0]), 'g-')
    pyplot.plot(x, y, ".k")
    pyplot.plot(x1, gauss(x1, coeff_first[0], coeff_first[1], coeff_first[2]), "-r")
    pyplot.fill_between(x1, ymin, gauss(x1, coeff_first[0], coeff_first[1], coeff_first[2]), color="r", alpha=0.1,
                        zorder=100)
    pyplot.plot(x2, gauss(x2, coeff_second[0], coeff_second[1], coeff_second[2]), "-b")
    pyplot.fill_between(x2, ymin, gauss(x2, coeff_second[0], coeff_second[1], coeff_second[2]), color="b", alpha=0.1,
                        zorder=101)
    pyplot.text(coeff_first[1], coeff_first[0] + 0.05 * yspan, "R2 = {0:.4f}".format(coeff_first[1]),
                horizontalalignment='center', verticalalignment='center')
    pyplot.text(coeff_second[1], coeff_second[0] + 0.05 * yspan, "R1 = {0:.4f}".format(coeff_second[1]),
                horizontalalignment='center', verticalalignment='center')
    pyplot.savefig("pruby.png", bbox_inches='tight', format="png")
    pyplot.clf()
    return


def ragan_r_line_estimation(temperature, peak_number):
    """	Input:	temperature in K (float)
                peak_number 1 for R1, 2 for R2 (int)
        Return:	wavenumber of ruby peak in the given temperature and room-pressure (float)
        Descr.: Calculate the value described above	"""
    t = temperature
    if peak_number == 1:
        return 14423. + 4.49 * 10 ** (-2) * t - 4.81 * 10 ** (-4) * t ** 2 + 3.71 * 10 ** (-7) * t ** 3
    elif peak_number == 2:
        return 14452. + 3.00 * 10 ** (-2) * t - 3.88 * 10 ** (-4) * t ** 2 + 2.55 * 10 ** (-7) * t ** 3
    else:
        return float("-Inf")


def seek_maxima(dots, method="primitive", smoothing_box_parameters=None):
    """	Input:	experimental x, y spectrum data (numpy n x 2 array)
                method of maxima seeking (str, "primitive" or "gaussian")
                [width of smoothing convolution box (int), convolution repetitions (int)]
        Return:	[[R2 peak pos. (float), R2 peak intensity (float),
                R2 error peak pos. (float), R2 error peak intensity (float)],
                [[R1 peak pos. (float), R1 peak intensity (float),
                R1 error peak pos. (float), R1 error peak intensity (float)]]
        Descr.: Find peaks, find their middles, order drawing them.	"""
    if smoothing_box_parameters == None:
        smoothing_box_parameters = [5, 5]
    maxima = []
    mean_intensity = np.mean(dots[:1])
    if method == "primitive":
        smoothdots = smoothen(dots, smoothing_box_parameters)
        for i in range(2, int(smoothdots.size / 2 - 2)):
            if (smoothdots[i, 1] >= np.max(smoothdots[i - 2:i + 2, 1] - 0.1) and
                        smoothdots[i, 0] > 692. and smoothdots[i, 0] < 700. and
                        smoothdots[i, 1] > mean_intensity):
                maxima.append([smoothdots[i, 0], smoothdots[i, 1],
                               abs(smoothdots[i + 1, 0] - smoothdots[i - 1, 0]) / 2., 0])
        maxima.append([float("-Inf"), float("-Inf"), float("-Inf"), float("-Inf")])
        maxima.append([float("-Inf"), float("-Inf"), float("-Inf"), float("-Inf")])
        return [maxima[0], maxima[1]]

    if method == "gaussian":
        maxima = seek_maxima(dots, method="primitive")
        ran = 0.20
        if (maxima[0][0] < 0 or maxima[1][0] < 0):
            return [maxima[0], maxima[1]]
        dots_first = np.empty(shape=[0, 2])
        dots_second = np.empty(shape=[0, 2])
        dots_for_drawing = np.empty(shape=[0, 2])
        for i in range(dots.shape[0]):
            if abs(dots[i, 0] - maxima[0][0]) < ran:
                dots_first = np.append(dots_first, [dots[i, 0:2]], axis=0)
            if abs(dots[i, 0] - maxima[1][0]) < ran:
                dots_second = np.append(dots_second, [dots[i, 0:2]], axis=0)
            if dots[i, 0] - maxima[0][0] > - 1. - ran and dots[i, 0] - maxima[1][0] < 1 + ran:
                dots_for_drawing = np.append(dots_for_drawing, [dots[i, 0:2]], axis=0)
        initial_guess_first = maxima[0][1], maxima[0][0], 0.35
        initial_guess_second = maxima[1][1], maxima[1][0], 0.35
        coeff_first, var_matrix_first = curve_fit(gauss, dots_first[:, 0], dots_first[:, 1], p0=initial_guess_first)
        coeff_second, var_matrix_second = curve_fit(gauss, dots_second[:, 0], dots_second[:, 1],
                                                    p0=initial_guess_second)
        coeff_err_first = np.sqrt(np.diag(var_matrix_first))
        coeff_err_second = np.sqrt(np.diag(var_matrix_second))
        maxima[0] = [coeff_first[1], coeff_first[0], coeff_err_first[1], coeff_err_first[0]]
        maxima[1] = [coeff_second[1], coeff_second[0], coeff_err_second[1], coeff_err_second[0]]
        plot_gaussian(dots_for_drawing, [maxima[0][0], maxima[1][0]], coeff_first, coeff_second, ran)
        return [maxima[0], maxima[1]]
    else:
        return [[float("-Inf"), float("-Inf"), float("-Inf"), float("-Inf")],
                [float("-Inf"), float("-Inf"), float("-Inf"), float("-Inf")]]


def smoothen(dots, smoothing_box_parameters):
    """	Input:	experimental x, y spectrum data (numpy n x 2 array)
                [width of smoothing convolution box (int), convolution repetitions (int)]
        Return:	smoothened experimental data (numpy n x 2 array)
        Descr.: Smoothen the data for further peak search	"""
    box_width = smoothing_box_parameters[0]
    box_times = smoothing_box_parameters[1]
    box = np.ones(box_width) / box_width
    for i in range(box_times):
        dots[:, 1] = np.convolve(dots[:, 1], box, mode='same')
    return dots


def wavenumber_to_wavelength(wavenumber):
    """	Input:	wavenumber in cm^-1 (float)
        Return:	wavelength in nm (float)
        Descr.: Calculate the value described above	"""
    return 10000000. / float(wavenumber)


class LogHandler:
    """ This thing somehow handles output and writing the log. """
    def __init__(self):
        self.time = datetime.datetime.now().strftime("%y-%m-%d")
        self.logpath = "pruby_{0}.log".format(self.time)
        self.logfile = open(self.logpath, "a")
        self.logfile.write(
            "APPENDING TO THE LOG FILE\nTIME: {0}\n".format
            (datetime.datetime.now().strftime("%Y/%m/%d %H:%M")))

    def reinitialise(self):
        self.logpath = "pruby_{0}.log".format(self.time)
        self.logfile = open(self.logpath, "a")
        self.logfile.write(
            "APPENDING TO THE LOG FILE\nTIME: {0}\n".format
            (datetime.datetime.now().strftime("%Y/%m/%d %H:%M")))

    def add(self, content):
        """	Input:	content (str)
            Return:	None
            Descr.:	Add content to the log	"""
        self.logfile.write(content + "\n")

    def close(self):
        """	Input:	content None
            Return:	None
            Descr.:	Close the log	"""
        temp = datetime.datetime.now().strftime("%Y/%m/%d %H:%M")
        self.logfile.write("CLOSING THE LOG FILE\nTIME: {}\n\n".format(temp))
        self.logfile.close()

    def ask(self, question=""):
        """	Input:	question (str)
            Return:	user input (str)
            Descr.:	Ask user about input in console, document in log.	"""
        temp = " < {0}".format(question)
        x = input(temp)
        self.add("{0}{1}".format(temp, str(x)))
        return x

    def say(self, text=""):
        """	Input:	text (str)
            Return:	None
            Descr.:	Write out the text to console and log	"""
        temp = " | {0}".format(str(text))
        print(temp)
        self.add(temp)
        return


log = LogHandler()


class PressureFromRuby:
    """	Descr.: A program for a God-knows-why reason placed inside an object. """
    def __init__(self):
        self.module = " "
        self.reference = [692.62, 694.03] # R2 and R1 values from one random ruby
        self.temperature = 293.15
        self.temperature_reference = 293.15
        self.temperature_sample = 293.15
        self.i()

    def menu(self):
        """	Descr.: Dish out a menu with module selection """
        self.module = log.ask("Choose a module: ")
        self.module = self.module + "h"
        self.module = self.module[0]
        if self.module in "rRpPtowhaq":
            getattr(self, self.module, None)()
        else:
            log.say("Module doesn't exist. Opening help:")
            self.h()
        return

    def a(self):
        """	Descr.: Write out the contributors. """
        log.say("PRUBY: Pressure from RUBY by Daniel Tchon")
        log.say("Version 0.1.170503")
        log.say("Pressure estimation based on Piermarini et al. (1975)")
        log.say("Temperature correction based on Ragan et al. (1992)")
        return

    def h(self):
        """	Descr.: Write out list of available modules """
        log.say("\tr - Reference change based on data file")
        log.say("\tR - Reference change based on manual input")
        log.say("\tp - Pressure estimation based on data file")
        log.say("\tP - Pressure estimation based on manual input")
        log.say("\tt - Manual change of current temperature")
        log.say("\tw - Change working directory ")
        log.say("\to - Overview of pruby temperature variables")
        log.say("\th - Help with module names")
        log.say("\ta - About the program")
        log.say("\tq - Quit the program")
        return

    def i(self):
        """	Descr.: Initiate the program (welcome + help) """
        log.say("PRUBY: Pressure from RUBY")
        self.h()
        return

    def o(self):
        """	Descr.: Write out local variables (mainly temperatures) """
        for var in self.__dict__:
            log.say("{0}: {1}".format(var, self.__dict__[var]))

    def p(self):
        """	Descr.: Calculate pressure based on spectrum. """
        self.temperature_sample = self.temperature
        path = log.ask("Input ruby file path: ")
        try:
            np.loadtxt(path)
        except Exception:
            log.say("File {} not found. Returning to menu.".format(path))
            return
        dots = np.loadtxt(path)
        maxima = seek_maxima(dots, "gaussian")
        if maxima[0][0] < 0. or maxima[1][0] < 0.:
            log.say("Not enough maxima found in the data.")
            log.say("Consider revising dataset or changing method.")
            return
        log.say("R1 =\t{0:.4f} +/- {1:.4f} nm".format(maxima[1][0], maxima[1][2]))
        log.say("R2 =\t{0:.4f} +/- {1:.4f} nm".format(maxima[0][0],maxima[0][2]))
        shift = [maxima[0][0] - self.reference[0], maxima[1][0] - self.reference[1]]
        pressure = [calculate_pressure_shift(shift[0], 2, self.temperature_reference, self.temperature_sample),
                    calculate_pressure_shift(shift[1], 1, self.temperature_reference, self.temperature_sample)]
        log.say("For T_ref = {0} and T_sample = {1}:".format(self.temperature_reference, self.temperature_sample))
        log.say("P1 =\t{0:.4f} GPa".format(pressure[1]))
        log.say("P2 =\t{0:.4f} GPa".format(pressure[0]))
        log.say("Pm =\t{0:.4f} GPa".format(pressure[0] / 2. + pressure[1] / 2.))
        return

    def P(self):
        """	Descr.: Calculate pressure based on manual input. """
        self.temperature_sample = self.temperature
        temp = []
        temp.append(float(log.ask("Input R1 peak position: ")))
        temp.append(float(log.ask("Input R2 peak position: ")))
        shift = [temp[1] - self.reference[0], temp[0] - self.reference[1]]
        pressure = [calculate_pressure_shift(shift[0], 2, self.temperature_reference ,self.temperature_sample),
                    calculate_pressure_shift(shift[1], 1, self.temperature_reference ,self.temperature_sample)]
        log.say("For T_ref = {0} and T_sample = {1}:".format(self.temperature_reference, self.temperature_sample))
        log.say("P1 =\t{0:.4f} GPa".format(pressure[1]))
        log.say("P2 =\t{0:.4f} GPa".format(pressure[0]))
        log.say("Pm =\t{0:.4f} GPa".format(pressure[0] / 2. + pressure[1] / 2.))
        return

    def q(self):
        """	Descr.: Save yourself and qut the program. """
        log.say("Program terminated. Log file saved as " + str(log.logpath))
        log.close()
        quit()

    def r(self):
        """	Descr.: Find reference peaks positions based on spectrum """
        self.temperature_reference = self.temperature
        path = log.ask("Input reference file path: ")
        try:
            np.loadtxt(path)
        except Exception:
            log.say("File {} not found. Returning to menu.".format(path))
            return
        dots = np.loadtxt(path)
        maxima = seek_maxima(dots, "gaussian")[-2:]
        if (maxima[0][0] < 0 or maxima[1][0] < 0):
            log.say("Not enough maxima found in the data.")
            log.say("Consider revising dataset or changing method.")
            return
        log.say("R1 =\t{0:.4f} +/- {1:.4f} nm".format(maxima[1][0], maxima[1][2]))
        log.say("R2 =\t{0:.4f} +/- {1:.4f} nm".format(maxima[0][0], maxima[0][2]))
        self.reference = [maxima[0][0], maxima[1][0]]
        return

    def R(self):
        """	Descr.: Find reference peaks positions based on manual input """
        self.temperature_reference = self.temperature
        temp = []
        temp.append(float(log.ask("Input R1 peak position: ")))
        temp.append(float(log.ask("Input R2 peak position: ")))
        log.say("R1 =\t{0:.4f} nm".format(temp[0]))
        log.say("R2 =\t{0:.4f} nm".format(temp[1]))
        self.reference = [temp[1],temp[0]]
        return

    def t(self):
        """	Descr.: Input current temperature (before p/P/r/R modules) """
        temp = float(log.ask("Input current temperature (in Celcius degrees): ")) + 273.15
        log.say("Temperature recalculated into Kelvins: " + "{0:.2f}".format(temp))
        self.temperature = temp
        return

    def w(self):
        """ Descr.: Change working directory for input and output (but not log)"""
        wd = log.ask("Specify new working directory: ")
        if os.path.isdir(wd):
            oldwd = os.getcwd()
            os.chdir(wd)
            log.say("Directory changed to {}.".format(os.getcwd()))
            log.close()
            log.reinitialise()
            log.add("LOG REINITIALISED FROM {}.".format(oldwd))
        else:
            log.say("Directory {} not found. Returning to menu.".format(wd))


# an object containing the program (lol)
pRuby = PressureFromRuby()
while True:
    pRuby.menu()
