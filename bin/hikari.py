import random
import numpy as np
import numpy.linalg as lin
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import pyqtplot for 3D - much faster than MatPlotLib


class Crystal:
    """Simple class for storing crystal parameters, content and properties.

    Space  | Direct [A] | Reciprocal [A] |
    Scalar | a_d        | a_r            |
    Vector | a_v        | a_w            |"""

    def __init__(self):
        self.a_d = 10.
        self.b_d = 15.
        self.c_d = 20.
        self.al_d = np.pi / 2
        self.be_d = np.pi / 2
        self.ga_d = np.pi / 2
        self.patterns = []

    def edit_cell(self, **parameters):
        """Edit unit cell using dictionary in Angstrom and degrees"""
        for key, value in parameters.items():
            key += '_d'
            if key in ('a_d', 'b_d', 'c_d'):
                setattr(self, key, value)
            elif key in ('alpha_d', 'beta_d', 'gamma_d'):
                setattr(self, key, np.radians(value))

    @property
    def v_d(self):
        """Provide a volume of unit cell in direct cell in Angstrom^3"""
        return self.a_d * self.b_d * self.c_d * \
               (1 - np.cos(self.al_d) ** 2 - np.cos(self.be_d) ** 2
                - np.cos(self.ga_d) ** 2 + 2 * np.cos(self.al_d)
                * np.cos(self.be_d) * np.cos(self.ga_d)) ** 0.5

    @property
    def a_v(self):
        """Provide an "a" vector in cartesian system in Angstrom"""
        return np.array((self.a_d, 0, 0))

    @property
    def b_v(self):
        """Provide a "b" vector in cartesian system in Angstrom"""
        return np.array((self.b_d * np.cos(self.ga_d),
                         self.b_d * np.sin(self.ga_d), 0))

    @property
    def c_v(self):
        """Provide a "c" vector in cartesian system in Angstrom"""
        return np.array((self.c_d * np.cos(self.be_d),
                         self.c_d * (np.cos(self.al_d) - np.cos(self.be_d)
                                     * np.cos(self.ga_d)) / np.sin(self.ga_d),
                         self.v_d / (self.a_d * self.b_d * np.sin(self.ga_d))))

    @property
    def a_r(self):
        """Provide an "a*" parameter in Angstrom^-1"""
        return 1 / self.a_d

    @property
    def b_r(self):
        """Provide a "b*" parameter in Angstrom^-1"""
        return 1 / self.b_d

    @property
    def c_r(self):
        """Provide a "c*" parameter in Angstrom^-1"""
        return 1 / self.c_d

    @property
    def al_r(self):
        """Provide an "alpha*" parameter in Angstrom^-1"""
        return np.pi - self.al_d

    @property
    def be_r(self):
        """Provide a "beta*" parameter in Angstrom^-1"""
        return np.pi - self.be_d

    @property
    def ga_r(self):
        """Provide a "gamma*" parameter in Angstrom^-1"""
        return np.pi - self.ga_d

    @property
    def v_r(self):
        """Provide a "v*" parameter in Angstrom^-3"""
        return 1 / self.v_d

    @property
    def a_w(self):
        """Provide an "a*" vector in reciprocal system in A"""
        return np.array((self.a_r, 0, 0))

    @property
    def b_w(self):
        """Provide a "b*" vector in reciprocal system in A"""
        return np.array((self.b_r * np.cos(self.ga_r),
                         self.b_r * np.sin(self.ga_r), 0))

    @property
    def c_w(self):
        """Provide a "c*" vector in reciprocal system in A"""
        return np.array((self.c_r * np.cos(self.be_r),
                      self.c_r * (np.cos(self.al_r) - np.cos(self.be_r)
                                  * np.cos(self.ga_r)) / np.sin(self.ga_r),
                      self.v_r / (self.a_r * self.b_r * np.sin(self.ga_r))))


class Pattern:
    """Single crystal diffraction pattern container"""
    reading_parameters = {'h': 4, 'k': 4, 'l': 4, 'i': 8, 's': 8, 'b': 4}
    writing_parameters = {'h': 4, 'k': 4, 'l': 4, 'i': 8, 's': 8, 'b': 4}
    writing_separator = False

    def __init__(self):
        self.wavelength = 0.71069
        self.reflections = []

    # keys: h, k, l - reciprocal lattice indexes
    # i - intensity/structure_factor, s - sigma,
    # b - batch/run, p - position

    def __add__(self, other):
        max_b = max(self.reflections, key=lambda x: x['b'])['b']
        sum_of_reflections = self.reflections
        for reflection in other.reflections:
            temp_reflection = reflection
            temp_reflection['b'] += max_b
            sum_of_reflections.append(temp_reflection)
        return sum_of_reflections

    @classmethod
    def writing_format(cls):
        """Format for writing the hkl file lines"""
        form = str('')
        for key in cls.writing_parameters.keys():
            form += '{' + key + ':>' + str(cls.writing_parameters[key]) + \
                    '.' + str(cls.writing_parameters[key] -
                              int(cls.writing_separator)) + '}'
        return form

    def edit_wavelength(self, wavelength):
        """Define wavelength [Ag, Mo, Cu or custom in A] used for measurement"""
        sources = {'Cr': 2.2896,	'Fe': 1.9360,	'Co': 1.79,
                   'Cu': 1.54056,	'Mo': 0.71069,	'Zr': 0.69,
                   'Ag': 0.5609}
        try:
            self.wavelength = sources[wavelength[:2]]
        except KeyError:
            self.wavelength = float(wavelength)
        except ValueError:
            pass

    def add_reflection(self, reflection):
        """Append specified reflection to the end of reflections list"""
        self.reflections.append(reflection)

    def delete_reflection(self, index):
        """Delete reflection with given index if exists"""
        try:
            del (self.reflections[index])
        except IndexError:
            pass

    def edit_reflection(self, index, **reflection):
        """Modify selected attributes of a reflection"""
        for key, value in reflection.items():
            self.reflections[index][key] = value

    def format_reflection(self, index):
        """Format reflection with given index using writing_format"""
        reflection_str = {}
        for key in Pattern.writing_parameters.keys():
            reflection_str[key] = str(self.reflections[index][key])
        return Pattern.writing_format().format(**reflection_str)

    def seek_reflection(self, **param):
        """Find indexes of all reflections described by dictionary of conditions"""
        positions = []
        for index, reflection in enumerate(self.reflections):
            if all(reflection[key] == param[key] for key in param.keys()):
                positions.append(index)
        return positions

    def clear_zero(self):
        reflections_to_delete = []
        for index, reflection in enumerate(self.reflections):
            if reflection['h'] == reflection['k'] == reflection['l'] == 0:
                reflections_to_delete.append(index)
        for index in reversed(reflections_to_delete):
            self.delete_reflection(index)

    def place(self, crystal):
        """Assign reflections their positions and delete h=k=l=0 reflections"""
        for index, reflection in enumerate(self.reflections):
            reflection['p'] = reflection['h'] * crystal.a_w + \
                              reflection['k'] * crystal.b_w + \
                              reflection['l'] * crystal.c_w

    def read(self, path):
        """Read .hkl file as specified by path and read_format"""
        hkl = open(path, 'r')
        for line in hkl:
            reflection = {'h': 0, 'k': 0, 'l': 0, 'i': 1., 's': 1., 'b': 1}
            for key in Pattern.reading_parameters.keys():
                try:
                    par = line[:Pattern.reading_parameters[key]]
                    if key in ('h', 'k', 'l', 'b'):
                        par = int(par)
                    elif key in ('i', 's'):
                        par = float(par)
                except ValueError:
                    continue
                reflection[key] = par
                line = line[Pattern.reading_parameters[key]:]
            self.add_reflection(reflection)
        hkl.close()

    def write(self, path):
        """Write .hkl file as specified by path and write_format"""
        hkl = open(path, 'w')
        for index, reflection in enumerate(self.reflections):
            hkl.write(self.format_reflection(index) + '\n')
        hkl.close()

    def dac(self, opening_angle=40, crystal_position=(1, 0, 0)):
        """Cut the reflections based on DAC angle (in angles)
        and [h, k, l] indexes of diamond-parallel crystal face"""
        opening_angle = np.radians([opening_angle])[0]
        normal = np.array(crystal_position)/lin.norm(np.array(crystal_position))
        reflections_to_delete = []
        for index, reflection in enumerate(self.reflections):
            proj = reflection['p'] - normal * np.dot(reflection['p'], normal)
            if np.isclose([lin.norm(proj)], [0.]):
                reflections_to_delete.append(index)
                continue
            phi = np.arccos(min(np.dot(proj, reflection['p']) / (lin.norm(proj)
                                * lin.norm(reflection['p'])), 1.))
            r_max = 2 * np.sin(max(opening_angle-phi, 0)) / self.wavelength
            if lin.norm(reflection['p']) > r_max:
                reflections_to_delete.append(index)
        for index in reversed(reflections_to_delete):
            self.delete_reflection(index)

    def thin_out(self, falloff=0.2):
        """Cut reflections based on a exp(-falloff*r) function"""
        reflections_to_delete = []
        for index, reflection in enumerate(self.reflections):
            r = (lin.norm(reflection['p']) ** 2)
            if random.random() > np.exp(-r * falloff):
                reflections_to_delete.append(index)
        for index in reversed(reflections_to_delete):
            self.delete_reflection(index)

    def trim(self, limit):
        """Cut the reflection further then the limit in A-1"""
        reflections_to_delete = []
        for index, reflection in enumerate(self.reflections):
            distance = lin.norm(reflection['p'])
            reflections_to_delete.append(index) if distance > limit else 0
        for index in reversed(reflections_to_delete):
            self.delete_reflection(index)

    @staticmethod
    def average(reflections):
        h, k, l = reflections[0]['h'], reflections[0]['k'], reflections[0]['l']
        b = len(reflections)
        i = sum(reflection['i'] for reflection in reflections) / b
        s = np.sqrt(
            sum(reflection['s'] ** 2 for reflection in reflections)) / b
        try:
            p = reflections[0]['p']
            return {'h': h, 'k': k, 'l': l, 'i': i, 's': s, 'b': b, 'p': p}
        except KeyError:
            return {'h': h, 'k': k, 'l': l, 'i': i, 's': s, 'b': b}

    def sort(self, key='i', descending=False):
        """Perform stable sort of reflections by a key"""
        self.reflections = sorted(self.reflections, key=lambda k: k[key],
                             reverse=descending)

    def reduce(self):
        """Average redundant reflections and write redundancy to batch number"""
        self.sort(key='l')
        self.sort(key='k')
        self.sort(key='h')
        same = [self.reflections.pop()]
        new_reflections = []
        while True:
            try:
                reflection = self.reflections.pop()
            except IndexError:
                new_reflections.append(self.average(same))
                break
            if all(reflection[key] == same[0][key] for key in ('h', 'k', 'l')):
                same.append(reflection)
            else:
                new_reflections.append(self.average(same))
                same = [reflection]
        self.reflections = new_reflections

    def draw(self, color_scheme='gist_rainbow', dpi=600, itosigma=False,
             legend=True, projection=('h', 'k', 0), property='batch', restrict=None,
             savefig=False, savename='hikari', scale=1.0, showfig=False):
        """Draw a cross-section of reciprocal lattice for given pattern

            color_scheme (colormap) Any colormap supported by matplotlib.cm
            dpi          (integer)  Dots Per Inch, quality of saved graphics
            itosigma     (boolean)  I/sigma(I) to be printed as transparency
            legend       (boolean)  Legend of used colors should be printed
            projection   (tuple)    cross-section ('h', 2, 'l') or None for 3D
            property     (string)   Colored property: [batch]/redundancy/none
            restrict     (tuple)    Restrict value of properties to certain set
            savefig      (boolean)  Save a figure to the file
            savename     (string)   Name of the saved figure
            scale        (float)    Scale factor for the reflection size
            showfig      (boolean)  Show a figure in matplotlib window"""

        # setting general geometry and parameters
        distance, axes, mode, direction, title = 1, [], '2d', 'h', 'hikari'
        backup_reflections = self.reflections
        fig = plt.figure()
        try:
            for value in projection:
                try:
                    distance = int(value)
                except ValueError:
                    axes.append(value)
        except TypeError:
                mode = '3d'
        if mode == '2d':
            direction = ({'h', 'k', 'l'} - set(axes)).pop()
            title = str(projection).replace(', ', '').replace('\'', '')
            ax = fig.add_subplot(1, 1, 1)
            axes.append('l')
        else:
            axes = ['h', 'k', 'l']
            ax = fig.add_subplot(111, projection='3d')

        # deciding about colored property
        if property == 'none':
            self.reduce()
            self.reduce()
        elif property == 'redundancy':
            self.reduce()

        # preparing for creating the lists of parameters
        x_pos, y_pos, z_pos, size, color, edge = [], [], [], [], [], []
        i_max = max(self.reflections, key=lambda x: x['i'])['i']
        itos_max = max(self.reflections, key=lambda x: x['i'] / x['s'])['i'] / \
                   max(self.reflections, key=lambda x: x['i'] / x['s'])['s']
        b_min = min(self.reflections, key=lambda x: x['b'])['b']
        b_max = max(self.reflections, key=lambda x: x['b'])['b']

        # preparing the colors
        b_range = b_max - b_min + 1
        color_map = cm.get_cmap(color_scheme)
        colors, alp = [], (1,)
        [colors.append(color_map(i / b_range)) for i in range(b_range)]

        # iterating over reflections
        for index, reflection in enumerate(self.reflections):
            if reflection[direction] is not distance and mode == '2d':
                continue
            if restrict is not None:
                if reflection['b'] not in restrict:
                    continue
            coordinates = {'h': 0, 'k': 1, 'l': 2}
            x_pos.append(reflection['p'][coordinates[axes[0]]])
            y_pos.append(reflection['p'][coordinates[axes[1]]])
            z_pos.append(reflection['p'][coordinates[axes[2]]])
            size.append(scale ** 2 * 100 * abs(reflection['i'] / i_max) ** 0.25)
            if itosigma:
                alp = abs(reflection['i'] / (reflection['s'] * itos_max))
                alp = (alp ** 0.25,)
            color.append(colors[(reflection['b'] - b_min)][:3] + alp)
            edge.append('None') if reflection['i'] > 0 else edge.append('k')

        # painting the plot
        directions = {'h': 'a* [A^-1]', 'k': 'b* [A^-1]', 'l': 'c* [A^-1]'}
        if mode == '2d':
            ax.set_title(title)
            ax.set_xlabel(directions[axes[0]])
            ax.set_ylabel(directions[axes[1]])
            ax.scatter(0, 0, s=20, c='k', marker='x')
            ax.scatter(x_pos, y_pos, s=size, c=color, marker='.',
                    edgecolors=edge, linewidth=[0.05 * s for s in size])
        elif mode == '3d':
            ax.scatter(0, 0, 0, s=20, c='k', marker='x')
            ax.scatter(x_pos, y_pos, z_pos, zdir='z', s=size,
                       c=color, marker='.')
            ax.set_xlabel(directions[axes[0]])
            ax.set_ylabel(directions[axes[1]])
            ax.set_zlabel(directions[axes[2]])
        ax.axis('equal')

        # preparing legend
        if legend and not property == 'none' and not b_range == 1:
            used_colors, used_batches = [], []
            color_skip = int(b_range / 25) + 1
            for i in range(0, b_range, color_skip):
                used_colors.append(plt.Rectangle((0, 0), 1, 1, fc=colors[i]))
                used_batches.append(str(i + b_min))
            if b_range < 10:
                ax.legend(used_colors, used_batches, loc=1,
                          prop={'size': 7}, title=property)
            elif b_range < 25:
                ax.legend(used_colors, used_batches, loc=1, ncol=2,
                          prop={'size': 6}, title=property)
            else:
                ax.legend(used_colors, used_batches, loc=1, ncol=3,
                          prop={'size': 5}, title=property)

        # returning reflections to their original state and saving figure
        self.reflections = backup_reflections
        fig = plt.gcf()
        if savefig:
            fig.savefig('{}.png'.format(savename), bbox_inches=None, dpi=dpi)
        if showfig:
            plt.show()


if __name__ == '__main__':
    p = Pattern()
    dox = Crystal()
    dox.edit_cell(a=11, b=12, c=16)
    p.read('small.hkl')
    p.place(dox)
    p.dac(opening_angle=40, crystal_position=(0, 0, 1))
    p.write('output.hkl')
    f = p.draw2d(projection=None, property='batch', scale=1.5, savefig=True, savename='hk0')
    plt.show()

# crystal = Crystal()
# p = Pattern()
# p.add_reflection({'h' : 100, 'k' : 100, 'l' : 100})
# p.place(crystal)
# print(p.reflections[0])
# crystal.edit_cell(a=10000, b=20000, c=30000)
# p.place(crystal)
# print(p.reflections[0])/home/dtchon

# TODO 3D call visualise and to pyqtplot
