from hikari.utility import angle2rad
import numpy as np
import numpy.linalg as lin


class BaseFrame:
    """
    This class stores and manipulates basic information present
    in majority of crystallographic information files such as unit cell
    parameters stored in scalars and vectors.

    BaseFrame utilises the following notation for stored attributes:

    - The name begins from a unit cell property we are interested in:
        - "a", "b", "c" are used to describe
          lattice constants *a*, *b* and *c* as well as
          lattice vectors **a**, **b**, and **c**,
        - "al", "be", "ga" are used to describe lattice constants
          *alpha*, *beta* and *gamma*,
          as in angle between vectors **a**, **b** and **c**,
        - "v" is used to describe the unit cell volume,
          a mixed product between vectors **a**, **b** and **c**,
        - "x", "y", "z" are used to describe normalised unit cell parameters,
          also denoted as unit cell directions,
          as in vector **z** being a vector **c** of length equal one.

    - The unit cell parameter symbol is then followed by an underscore "_".

    - The name is ended by a single letter which describes
      if we are working in direct or reciprocal space,
      as well as if we want to access a scalar or a vector.

        - "d" (from Direct) is used to denote direct space scalars,
        - "r" (from Reciprocal) is used to denote reciprocal space scalars,
        - "v" (from Vector) is used to denote direct space vectors,
        - 'w" (similar to "v") is used to denote reciprocal space vectors.

    The values can be accessed by referencing a given attribute in the object,
    for example :class:`BaseFrame`. :attr:`a_d` stores information about
    the lattice constant *a* in direct space as a floating point,
    but :class:`BaseFrame`. :attr:`a_v` is a direct space vector. Please mind
    that only the meaningful combinations of descriptors are defined.

    Available attributes have been once again presented in a table below:

    +----------+------------+------------------+------------------+------------+
    |          | Available  | in direct        | in reciprocal    |Unit (^-1 in|
    |          | constants  | space            | space            |reciprocal) |
    +==========+============+==================+==================+============+
    | Scalars  | a, b, c    | a_d, b_d, c_d    | a_r, b_r, c_r    | Angstrom   |
    |          +------------+------------------+------------------+------------+
    |          | al, be, ga | al_d, be_d, ga_d | al_r, be_r, ga_r | Radian     |
    |          +------------+------------------+------------------+------------+
    |          | v          | v_d              | v_r              | Angstrom^3 |
    +----------+------------+------------------+------------------+------------+
    | Vectors  | a, b, c    | a_v, b_v, c_v    | a_w, b_w, c_w    | Angstrom   |
    |          +------------+------------------+------------------+------------+
    |          | x, y, z    | x_v, y_v, z_v    | x_w, y_w, z_w    | Angstrom   |
    +----------+------------+------------------+------------------+------------+
    """

    def __init__(self):
        self.__a_d = self.__b_d = self.__c_d = 1.0
        self.__a_r = self.__b_r = self.__c_r = 1.0
        self.__al_d = self.__be_d = self.__ga_d = np.pi / 2
        self.__al_r = self.__be_r = self.__ga_r = np.pi / 2
        self.__a_v, self.__b_v, self.__c_v = np.eye(3)
        self.__a_w, self.__b_w, self.__c_w = np.eye(3)
        self.edit_cell(a=1.0, b=1.0, c=1.0, al=90.0, be=90.0, ga=90.0)
        self.orientation = np.array(((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0)))
        """3x3 matrix describing orientation of crystal during experiment."""

    def edit_cell(self, **parameters):
        """
        Edit direct space unit cell using a dictionary.

        The input dictionary accepts only the following six keys:

        - "a" - for unit cell parameter *a* given in Angstrom,

        - "b" - for unit cell parameter *b* given in Angstrom,

        - "c" - for unit cell parameter *c* given in Angstrom,

        - "al" - for unit cell parameter *alpha* given in degrees or radians,

        - "be" - for unit cell parameter *beta* given in degrees or radians,

        - "ga" - for unit cell parameter *gamma* given in degrees or radians.

        This method is equivalent to manually setting all six unit cell
        parameters in direct space, :attr:`a_d`, :attr:`b_d`, :attr:`c_d`,
        :attr:`al_d`, :attr:`be_d`, :attr:`ga_d`, and then running a private
        method :meth:`_refresh_cell` to update other values.

        Please mind that the while the "a", "b" and "c" are always given in
        Angstrom, the angles might be given either in degrees or in radians.
        For further details, please consult
        :func:`hikari.utility.math_tools.both2rad` function.

        It is not required for all previously stated keys to be present
        at each method call. If a key has not been given, previously provided
        and stored value is being used. If no value has been given,
        the default length values of 1.0 for *a*, *b*, *c* and default angle
        values of *pi/2* for *al*, *be*, *ga* are used instead.

        :param parameters: Values of unit cell parameters to be changed
        :type parameters: float
        """
        for key, value in parameters.items():
            assert key in ('a', 'b', 'c', 'al', 'be', 'ga'), "unknown parameter"
            setattr(self, '{}_d'.format(key), value)
        self._refresh_cell()

    def from_cif_frame(self, frame):
        """
        Attempt importing unit cell parameters and orientation matrix
        from provided :class:`hikari.dataframes.cif.CifFrame` object.

        This method requires at least cell parameters to be defined
        in CifFrame object. It also attempts to import the orientation matrix,
        but passes if unsuccessful.

        :param frame: CifFrame containing cell parameters and,
            optionally, crystal orientation matrix.
        :type frame: kesshou.dataframes.CifFrame
        """

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

    def _refresh_cell(self):
        """
        Recalculate all vectors and scalars other than :attr:`a_d`,
        :attr:`b_d`, :attr:`c_d`, :attr:`al_d`, :attr:`be_d`, :attr:`ga_d`
        based on the currently stored values of the aforementioned six.
        """

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

    @property
    def a_d(self):
        """
        Scalar length of unit cell vector **a** in direct space.

        :return: Length of **a**.
        :rtype: float
        """
        return self.__a_d

    @a_d.setter
    def a_d(self, value):
        self.__a_d = value

    @property
    def b_d(self):
        """
        Scalar length of unit cell vector **b** in direct space.

        :return: Length of **b**.
        :rtype: float
        """
        return self.__b_d

    @b_d.setter
    def b_d(self, value):
        self.__b_d = value

    @property
    def c_d(self):
        """
        Scalar length of unit cell vector **c** in direct space.

        :return: Length of **c**.
        :rtype: float
        """
        return self.__c_d

    @c_d.setter
    def c_d(self, value):
        self.__c_d = value

    @property
    def al_d(self):
        """
        Scalar angle between unit cell vectors **b** and **c** in direct space.
        Setting this value will perform degrees / radian check according to
        documentation of :func:`hikari.utility.math_tools.both2rad`.

        :return: Angle between **b** and **c**.
        :rtype: float
        """
        return self.__al_d

    @al_d.setter
    def al_d(self, value):
        self.__al_d = angle2rad(value)

    @property
    def be_d(self):
        """
        Scalar angle between unit cell vectors **c** and **a** in direct space.
        Setting this value will perform degrees / radian check according to
        documentation of :func:`hikari.utility.math_tools.both2rad`.

        :return: Angle between **c** and **a**.
        :rtype: float
        """
        return self.__be_d

    @be_d.setter
    def be_d(self, value):
        self.__be_d = angle2rad(value)

    @property
    def ga_d(self):
        """
        Scalar angle between unit cell vectors **c** and **a** in direct space.
        Setting this value will perform degrees / radian check according to
        documentation of :func:`hikari.utility.math_tools.both2rad`.

        :return: Angle between **c** and **a**.
        :rtype: float
        """
        return self.__ga_d

    @ga_d.setter
    def ga_d(self, value):
        self.__ga_d = angle2rad(value)

    @property
    def v_d(self):
        """
        Mixed product of vectors **a**, **b** and **c** in direct space.

        :return: Unit cell volume in direct space.
        :rtype: float
        """
        return self.__v_d

    @property
    def a_v(self):
        """
        Unit cell vector **a** in direct space.

        :return: Unit cell vector **a** in direct space.
        :rtype: numpy.array
        """
        return self.__a_v

    @property
    def b_v(self):
        """
        Unit cell vector **b** in direct space.

        :return: Unit cell vector **b** in direct space.
        :rtype: numpy.array
        """
        return self.__b_v

    @property
    def c_v(self):
        """
        Unit cell vector **c** in direct space.

        :return: Unit cell vector **c** in direct space.
        :rtype: numpy.array
        """
        return self.__c_v

    @property
    def x_v(self):
        """
        Crystallographic direction vector **x** in direct space.

        :return: Normalised vector **a** in direct space.
        :rtype: numpy.array
        """
        return self.__a_v / lin.norm(self.__a_v)

    @property
    def y_v(self):
        """
        Crystallographic direction vector **y** in direct space.

        :return: Normalised vector **b** in direct space.
        :rtype: numpy.array
        """
        return self.__b_v / lin.norm(self.__b_v)

    @property
    def z_v(self):
        """
        Crystallographic direction vector **z** in direct space.

        :return: Normalised vector **c** in direct space.
        :rtype: numpy.array
        """
        return self.__c_v / lin.norm(self.__c_v)

    @property
    def a_r(self):
        """
        Scalar length of unit cell vector **a** in reciprocal space.

        :return: Length of **a\***.
        :rtype: float
        """
        return self.__a_r

    @a_r.setter
    def a_r(self, value):
        self.__a_r = value

    @property
    def b_r(self):
        """
        Scalar length of unit cell vector **b** in reciprocal space.

        :return: Length of **b\***.
        :rtype: float
        """
        return self.__b_r

    @b_r.setter
    def b_r(self, value):
        self.__b_r = value

    @property
    def c_r(self):
        """
        Scalar length of unit cell vector **c** in reciprocal space.

        :return: Length of **c\***.
        :rtype: float
        """
        return self.__c_r

    @c_r.setter
    def c_r(self, value):
        self.__c_r = value

    @property
    def al_r(self):
        """
        Scalar angle between unit cell vectors **b\*** and **c\***
        in reciprocal space.


        :return: Angle between **b\*** and **c\***.
        :rtype: float
        """
        return self.__al_r

    @al_r.setter
    def al_r(self, value):
        self.__al_r = angle2rad(value)

    @property
    def be_r(self):
        """
        Scalar angle between unit cell vectors **c\*** and **a\***
        in reciprocal space.

        :return: Angle between **c\*** and **a\***.
        :rtype: float
        """
        return self.__be_r

    @be_r.setter
    def be_r(self, value):
        self.__be_r = angle2rad(value)

    @property
    def ga_r(self):
        """
        Scalar angle between unit cell vectors **a\*** and **b\***
        in reciprocal space.

        :return: Angle between **a\*** and **b\***.
        :rtype: float
        """
        return self.__ga_r

    @ga_r.setter
    def ga_r(self, value):
        self.__ga_r = angle2rad(value)

    @property
    def v_r(self):
        """
        Mixed product of vectors **a\***, **b\*** and **c\***
        in reciprocal space.

        :return: Unit cell volume in direct space.
        :rtype: float
        """
        return self.__v_r

    @property
    def a_w(self):
        """
        Unit cell vector **a\*** in reciprocal space.

        :return: Unit cell vector **a\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self.__a_w

    @property
    def b_w(self):
        """
        Unit cell vector **b\*** in reciprocal space.

        :return: Unit cell vector **b\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self.__b_w

    @property
    def c_w(self):
        """
        Unit cell vector **c\*** in reciprocal space.

        :return: Unit cell vector **c\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self.__c_w

    @property
    def x_w(self):
        """
        Crystallographic direction vector **x\*** in reciprocal space.

        :return: Normalised vector **a\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self.__a_w / lin.norm(self.__a_w)

    @property
    def y_w(self):
        """
        Crystallographic direction vector **y\*** in reciprocal space.

        :return: Normalised vector **b\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self.__b_w / lin.norm(self.__b_w)

    @property
    def z_w(self):
        """
        Crystallographic direction vector **z\*** in reciprocal space.

        :return: Normalised vector **c\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self.__c_w / lin.norm(self.__c_w)
