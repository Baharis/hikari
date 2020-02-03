from kesshou.utility import angle
import numpy as np
import numpy.linalg as lin


class BaseFrame:
    """Base container for data common in Kesshou containers, eg.: unit cell

    Space          | Direct        | Reciprocal    |
    Scalar         | *_d           | *_r           |
    Vector         | *_v           | *_w           |

    a, b, c        | unit cell lengths in Angstrom / Angstrom^-1
    x, y, z        | normalised unit cell vectors
    al, be, ga     | unit cell angles in radians
    v              | unit cell volume in Angstrom^3 / Angstrom^-3"""

    def __init__(self):
        self.__a_d = self.__b_d = self.__c_d = 1.0
        self.__a_r = self.__b_r = self.__c_r = 1.0
        self.__al_d = self.__be_d = self.__ga_d = np.pi / 2
        self.__al_r = self.__be_r = self.__ga_r = np.pi / 2
        self.__a_v, self.__b_v, self.__c_v = np.eye(3)
        self.__a_w, self.__b_w, self.__c_w = np.eye(3)
        self.edit_cell(a=1.0, b=1.0, c=1.0, al=90.0, be=90.0, ga=90.0)
        self.orientation = np.array(((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0)))

    def edit_cell(self, **parameters):
        """Edit direct space unit cell using dictionary with the following keys:
        a, b, c [in Angstrom] and al, be, ga [in degrees or radians]."""

        # INSERT NEW VALUES OF DIRECT CELL PARAMETERS
        for key, value in parameters.items():
            assert key in ('a', 'b', 'c', 'al', 'be', 'ga'), "unknown parameter"
            setattr(self, '{}_d'.format(key), value)
        self._refresh_cell()

    def _refresh_cell(self):
        """A function to refresh all vectors and scalars of unit cell"""

        def calculate_volume(a, b, c, al, be, ga):
            return a * b * c * \
               (1 - np.cos(al) ** 2 - np.cos(be) ** 2 - np.cos(ga) ** 2
                + 2 * np.cos(al) * np.cos(be) * np.cos(ga)) ** 0.5
        self.__v_d = calculate_volume(self.a_d, self.b_d, self.c_d,
                                      self.al_d, self.be_d, self.ga_d)

        def calculate_reciprocal_cell(a, b, c, al, be, ga, v):
            a_r = b * c * np.sin(al) / v
            b_r = c * a * np.sin(be) / v
            c_r = a * b * np.sin(ga) / v
            al_r = np.arccos((np.cos(be) * np.cos(ga) - np.cos(al)) /
                             (np.sin(be) * np.sin(ga)))
            be_r = np.arccos((np.cos(ga) * np.cos(al) - np.cos(be)) /
                             (np.sin(ga) * np.sin(al)))
            ga_r = np.arccos((np.cos(al) * np.cos(be) - np.cos(ga)) /
                             (np.sin(al) * np.sin(be)))
            return a_r, b_r, c_r, al_r, be_r, ga_r
        self.a_r, self.b_r, self.c_r, self.al_r, self.be_r, self.ga_r = \
            calculate_reciprocal_cell(self.a_d, self.b_d, self.c_d,
                                      self.al_d, self.be_d, self.ga_d, self.v_d)
        self.__v_r = calculate_volume(self.a_r, self.b_r, self.c_r,
                                      self.al_r, self.be_r, self.ga_r)

        def calculate_vectors():
            self.__a_v = np.array((self.a_d, 0, 0))
            self.__b_v = np.array((self.b_d * np.cos(self.ga_d),
                                   self.b_d * np.sin(self.ga_d), 0))
            self.__c_v = np.array((self.c_d * np.cos(self.be_d),
                                   self.c_d * (
                                   np.cos(self.al_d) - np.cos(self.be_d)
                                   * np.cos(self.ga_d)) / np.sin(self.ga_d),
                                   self.v_d / (
                                   self.a_d * self.b_d * np.sin(self.ga_d))))
            self.__a_w = np.cross(self.b_v, self.c_v) / self.v_d
            self.__b_w = np.cross(self.c_v, self.a_v) / self.v_d
            self.__c_w = np.cross(self.a_v, self.b_v) / self.v_d
        calculate_vectors()

    def from_cif_frame(self, frame):
        """Import necessary crystal parameters from CifFrame"""

        # IMPORT AND CHANGE LATTICE PARAMETERS
        new_parameters = {
            'a': float(frame.data['_cell_length_a']),
            'b': float(frame.data['_cell_length_b']),
            'c': float(frame.data['_cell_length_c']),
            'al': float(frame.data['_cell_angle_alpha']),
            'be': float(frame.data['_cell_angle_beta']),
            'ga': float(frame.data['_cell_angle_gamma'])}
        self.edit_cell(**new_parameters)

        # IMPORT AND CHANGE ORIENTATION MATRIX
        try:
            self.orientation = \
                np.array(((float(frame.data['_diffrn_orient_matrix_UB_11']),
                           float(frame.data['_diffrn_orient_matrix_UB_12']),
                           float(frame.data['_diffrn_orient_matrix_UB_13'])),
                          (float(frame.data['_diffrn_orient_matrix_UB_21']),
                           float(frame.data['_diffrn_orient_matrix_UB_22']),
                           float(frame.data['_diffrn_orient_matrix_UB_23'])),
                          (float(frame.data['_diffrn_orient_matrix_UB_31']),
                           float(frame.data['_diffrn_orient_matrix_UB_32']),
                           float(frame.data['_diffrn_orient_matrix_UB_33']))))
        except KeyError:
            pass

    @property
    def a_d(self):
        return self.__a_d

    @a_d.setter
    def a_d(self, value):
        self.__a_d = value

    @property
    def b_d(self):
        return self.__b_d

    @b_d.setter
    def b_d(self, value):
        self.__b_d = value

    @property
    def c_d(self):
        return self.__c_d

    @c_d.setter
    def c_d(self, value):
        self.__c_d = value

    @property
    def al_d(self):
        return self.__al_d

    @al_d.setter
    def al_d(self, value):
        self.__al_d = angle(value)

    @property
    def be_d(self):
        return self.__be_d

    @be_d.setter
    def be_d(self, value):
        self.__be_d = angle(value)

    @property
    def ga_d(self):
        return self.__ga_d

    @ga_d.setter
    def ga_d(self, value):
        self.__ga_d = angle(value)

    @property
    def v_d(self):
        return self.__v_d

    @property
    def a_v(self):
        return self.__a_v

    @property
    def b_v(self):
        return self.__b_v

    @property
    def c_v(self):
        return self.__c_v

    @property
    def x_v(self):
        return self.__a_v / lin.norm(self.__a_v)

    @property
    def y_v(self):
        return self.__b_v / lin.norm(self.__b_v)

    @property
    def z_v(self):
        return self.__c_v / lin.norm(self.__c_v)

    @property
    def a_r(self):
        return self.__a_r

    @a_r.setter
    def a_r(self, value):
        self.__a_r = value

    @property
    def b_r(self):
        return self.__b_r

    @b_r.setter
    def b_r(self, value):
        self.__b_r = value

    @property
    def c_r(self):
        return self.__c_r

    @c_r.setter
    def c_r(self, value):
        self.__c_r = value

    @property
    def al_r(self):
        return self.__al_r

    @al_r.setter
    def al_r(self, value):
        self.__al_r = angle(value)

    @property
    def be_r(self):
        return self.__be_r

    @be_r.setter
    def be_r(self, value):
        self.__be_r = angle(value)

    @property
    def ga_r(self):
        return self.__ga_r

    @ga_r.setter
    def ga_r(self, value):
        self.__ga_r = angle(value)

    @property
    def v_r(self):
        return self.__v_r

    @property
    def a_w(self):
        return self.__a_w

    @property
    def b_w(self):
        return self.__b_w

    @property
    def c_w(self):
        return self.__c_w

    @property
    def x_w(self):
        return self.__a_w / lin.norm(self.__a_w)

    @property
    def y_w(self):
        return self.__b_w / lin.norm(self.__b_w)

    @property
    def z_w(self):
        return self.__c_w / lin.norm(self.__c_w)
