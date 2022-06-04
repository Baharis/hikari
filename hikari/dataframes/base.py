import numpy as np
from hikari.utility import angle2rad, cfloat, det3x3


class BaseFrame:
    """
    This class stores and manipulates basic information present
    in majority of crystallographic information files such as unit cell
    parameters stored in scalars and vectors.

    BaseFrame utilises the following notation for stored attributes:

    - The name begins from a unit cell property we are interested in:

     - "a", "b", "c" describe unit cell lengths/vectors *a*, *b*, *c*,
     - "al", "be", "ga" describe unit cell angles *alpha*, *beta*, *gamma*,
     - "v" describes unit cell volume,
     - "x", "y", "z" describe directions - normalised unit cell vectors.
     - "A", "G" describe stacked vector and metric matrix, respectively.

    - The unit cell parameter symbol is then followed by an underscore "_".

    - The name ends with a single letter denoting type of space and variable:

     - "d" (from Direct) denotes direct space scalars/matrices,
     - "r" (from Reciprocal) denotes reciprocal space scalars/matrices,
     - "v" (from Vector) denotes direct space vectors,
     - 'w" (similar to "v") denotes reciprocal space vectors.

    The values can be accessed by referencing a given attribute in the object,
    for example :class:`BaseFrame`. :attr:`a_d` stores information about
    the lattice constant *a* in direct space as a floating point,
    but :class:`BaseFrame`. :attr:`a_v` is a direct space vector.
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
    | Matrices | A          | A_d              | A_r              | Angstrom^2 |
    |          +------------+------------------+------------------+------------+
    |          | G          | G_d              | G_r              | Angstrom^2 |
    +----------+------------+------------------+------------------+------------+
    """

    IMPORTED_FROM_CIF = {
        'a': ['_cell_length_a', cfloat, 1.0],
        'b': ['_cell_length_b', cfloat, 1.0],
        'c': ['_cell_length_c', cfloat, 1.0],
        'al': ['_cell_length_alpha', cfloat, 90],
        'be': ['_cell_length_beta', cfloat, 90],
        'ga': ['_cell_length_gamma', cfloat, 90],
        'ub11': ['_diffrn_orient_matrix_UB_11', float, 1.0],
        'ub12': ['_diffrn_orient_matrix_UB_12', float, 0.0],
        'ub13': ['_diffrn_orient_matrix_UB_13', float, 0.0],
        'ub21': ['_diffrn_orient_matrix_UB_21', float, 0.0],
        'ub22': ['_diffrn_orient_matrix_UB_22', float, 1.0],
        'ub23': ['_diffrn_orient_matrix_UB_23', float, 0.0],
        'ub31': ['_diffrn_orient_matrix_UB_31', float, 0.0],
        'ub32': ['_diffrn_orient_matrix_UB_32', float, 0.0],
        'ub33': ['_diffrn_orient_matrix_UB_33', float, 1.0]}

    def __init__(self):
        self._a_d = self._b_d = self._c_d = 1.0
        self._a_r = self._b_r = self._c_r = 1.0
        self._al_d = self._be_d = self._ga_d = np.pi / 2
        self._al_r = self._be_r = self._ga_r = np.pi / 2
        self._a_v, self._b_v, self._c_v = np.eye(3)
        self._a_w, self._b_w, self._c_w = np.eye(3)
        self._refresh_cell()
        self.orientation = np.array(((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0)))
        """3x3 matrix describing orientation of crystal during experiment."""

    def edit_cell(self, **parameters):
        """
        Edit direct space unit cell using a dictionary with the following keys:

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
        For details see function :func:`hikari.utility.math_tools.angle2rad`.

        It is not required for all previously stated keys to be present
        at each method call. If a key has not been given, previously provided
        and stored value is being used. If no value has been given,
        the default length values of 1.0 for *a*, *b*, *c* and default angle
        values of *pi/2* for *al*, *be*, *ga* are used instead.

        :param parameters: Values of unit cell parameters to be changed
        :type parameters: float
        """
        for key, value in parameters.items():
            if key not in {'a', 'b', 'c', 'al', 'be', 'ga'}:
                raise KeyError(f'Unknown unit cell parameter: {key}')
            elif key in {'a', 'b', 'c'}:
                setattr(self, f'_{key}_d', value)
            elif key in {'al', 'be', 'ga'}:
                setattr(self, f'_{key}_d', angle2rad(value))
        self._refresh_cell()

    def fill_from_cif_block(self, block, fragile=False):
        """
        Import all data specified in :attr:`~.BaseFrame.IMPORTED_FROM_CIF`
        such as unit cell parameters and orientation matrix from provided
        instance of :class:`hikari.dataframes.cif.CifBlock` called `block`.
        Unless `fragile` is `True`, use defaults instead of rising KeyError.

        :param block: CifBlock containing imported information.
        :type block: hikari.dataframes.CifBlock
        :param fragile: If True, raise Error when any imported info is missing
        :type fragile: bool
        """
        imp = {}
        for k, v in self.IMPORTED_FROM_CIF.items():
            if fragile:
                imp[k] = block.get_as_type(v[0], v[1])
            else:
                imp[k] = block.get_as_type(v[0], v[1], v[2])
        cell_par_names = {'a', 'b', 'c', 'al', 'be', 'ga'}
        new_parameters = {k: v for k, v in imp.items() if k in cell_par_names}
        self.edit_cell(**new_parameters)
        orientation_value_names = ('ub11', 'ub12', 'ub13',
                                   'ub21', 'ub22', 'ub23',
                                   'ub31', 'ub32', 'ub33')
        orientation_values = [imp[n] for n in orientation_value_names]
        self.orientation = np.array(orientation_values).reshape(3, 3)

    def _refresh_cell(self):
        """
        Recalculate all vectors and scalars other than :attr:`a_d`,
        :attr:`b_d`, :attr:`c_d`, :attr:`al_d`, :attr:`be_d`, :attr:`ga_d`
        based on the currently stored values of the aforementioned six.
        """
        a, b, c = self.a_d, self.b_d, self.c_d
        sa, sb, sg = np.sin(self.al_d), np.sin(self.be_d), np.sin(self.ga_d)
        ca, cb, cg = np.cos(self.al_d), np.cos(self.be_d), np.cos(self.ga_d)
        v = a * b * c * np.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)

        self._a_v = np.array([a, 0, 0])
        self._b_v = np.array([b * cg, b * sg, 0])
        self._c_v = np.array([c * cb, c * (ca - cb * cg) / sg, v / (a * b * sg)])

        self._a_r = b * c * sa / v
        self._b_r = c * a * sb / v
        self._c_r = a * b * sg / v
        self._al_r = np.arccos((cb * cg - ca) / (sb * sg))
        self._be_r = np.arccos((cg * ca - cb) / (sg * sa))
        self._ga_r = np.arccos((ca * cb - cg) / (sa * sb))

        self._a_w = np.cross(self.b_v, self.c_v) / v
        self._b_w = np.cross(self.c_v, self.a_v) / v
        self._c_w = np.cross(self.a_v, self.b_v) / v

    @property
    def a_d(self):
        """
        :return: Length of unit cell vector **a** in direct space.
        :rtype: float
        """
        return self._a_d

    @property
    def b_d(self):
        """
        :return: Length of unit cell vector **b** in direct space.
        :rtype: float
        """
        return self._b_d

    @property
    def c_d(self):
        """
        :return: Length of unit cell vector **c** in direct space.
        :rtype: float
        """
        return self._c_d

    @property
    def al_d(self):
        """
        :return: Angle between vectors **b** and **c** in degrees.
        :rtype: float
        """
        return self._al_d

    @property
    def be_d(self):
        """
        :return: Angle between vectors **c** and **a** in degrees.
        :rtype: float
        """
        return self._be_d

    @property
    def ga_d(self):
        """
        :return: Angle between vectors **a** and **b** in degrees.
        :rtype: float
        """
        return self._ga_d

    @property
    def v_d(self):
        """
        :return: Unit cell volume in direct space.
        :rtype: float
        """
        return det3x3(self.A_d)

    @property
    def a_v(self):
        """
        :return: Unit cell vector **a** in direct space.
        :rtype: numpy.array
        """
        return self._a_v

    @property
    def b_v(self):
        """
        :return: Unit cell vector **b** in direct space.
        :rtype: numpy.array
        """
        return self._b_v

    @property
    def c_v(self):
        """
        :return: Unit cell vector **c** in direct space.
        :rtype: numpy.array
        """
        return self._c_v

    @property
    def A_d(self):
        """
        :return: Matrix A with vertically stacked direct space vectors.
        :rtype: np.array
        """
        return np.vstack([self._a_v, self._b_v, self._c_v])

    @property
    def G_d(self):
        """
        :return: Direct space metric matrix [ai . aj]ij.
        :rtype: np.array
        """
        return self.A_d @ self.A_d.T

    @property
    def a_r(self):
        """
        :return: Length of unit cell vector **a\*** in reciprocal space.
        :rtype: float
        """
        return self._a_r

    @property
    def b_r(self):
        """
        :return: Length of unit cell vector **b\*** in reciprocal space.
        :rtype: float
        """
        return self._b_r

    @property
    def c_r(self):
        """
        :return: Length of unit cell vector **c\*** in reciprocal space.
        :rtype: float
        """
        return self._c_r

    @property
    def al_r(self):
        """
        :return: Angle between vectors **b\*** and **c\*** in degrees.
        :rtype: float
        """
        return self._al_r

    @property
    def be_r(self):
        """
        :return: Angle between vectors **c\*** and **a\*** in degrees.
        :rtype: float
        """
        return self._be_r

    @property
    def ga_r(self):
        """
        :return: Angle between vectors **a\*** and **b\*** in degrees.
        :rtype: float
        """
        return self._ga_r

    @property
    def v_r(self):
        """
        :return: Unit cell volume in reciprocal space.
        :rtype: float
        """
        return det3x3(self.A_r)

    @property
    def a_w(self):
        """
        :return: Unit cell vector **a\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self._a_w

    @property
    def b_w(self):
        """
        :return: Unit cell vector **b\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self._b_w

    @property
    def c_w(self):
        """
        :return: Unit cell vector **c\*** in reciprocal space.
        :rtype: numpy.array
        """
        return self._c_w

    @property
    def A_r(self):
        """
        :return: Matrix A\* with vertically stacked reciprocal space vectors.
        :rtype: np.array
        """
        return np.vstack([self._a_w, self._b_w, self._c_w])

    @property
    def G_r(self):
        """
        :return: Reciprocal space metric matrix [ai\* . aj\*]ij.
        :rtype: np.array
        """
        return self.A_r @ self.A_r.T
