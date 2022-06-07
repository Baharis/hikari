import numpy as np
from uncertainties import ufloat, ufloat_fromstr, unumpy
from uncertainties.umath import sin as usin, cos as ucos
from uncertainties.umath import acos as uacos, sqrt as usqrt

from hikari.dataframes import BaseFrame


class UBaseFrame(BaseFrame):
    """
    A sub-class of :class:`hikari.dataframes.BaseFrame` capable of the same
    operation as its parent, but using `uncertainty.ufloats` instead of floats.
    As a result, types specified in docstring might be wrong due to inheritance.
    """

    IMPORTED_FROM_CIF = {
        'a': ['_cell_length_a', ufloat_fromstr, 1.0],
        'b': ['_cell_length_b', ufloat_fromstr, 1.0],
        'c': ['_cell_length_c', ufloat_fromstr, 1.0],
        'al': ['_cell_length_alpha', ufloat_fromstr, 90],
        'be': ['_cell_length_beta', ufloat_fromstr, 90],
        'ga': ['_cell_length_gamma', ufloat_fromstr, 90],
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
        super(UBaseFrame, self).__init__()
        u0, u1, upi = ufloat(0., 0), ufloat(1., 0), ufloat(np.pi, 0)
        u_eye = unumpy.uarray(np.eye(3), np.zeros([3, 3]))
        self._a_d = self._b_d = self._c_d = u1
        self._a_r = self._b_r = self._c_r = u1
        self._al_d = self._be_d = self._ga_d = upi / 2
        self._al_r = self._be_r = self._ga_r = upi / 2
        self._a_v, self._b_v, self._c_v = u_eye
        self._a_w, self._b_w, self._c_w = u_eye
        self._refresh_cell()
        self.orientation = np.array(((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0)))
        """3x3 matrix describing orientation of crystal during experiment."""

    def _refresh_cell(self):
        a, b, c = self.a_d, self.b_d, self.c_d
        sa, sb, sg = usin(self.al_d), usin(self.be_d), usin(self.ga_d)
        ca, cb, cg = ucos(self.al_d), ucos(self.be_d), ucos(self.ga_d)
        v = a * b * c * usqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg)

        u0 = ufloat(0., 0)
        self._a_v = np.array([a, u0, u0])
        self._b_v = np.array([b * cg, b * sg, u0])
        self._c_v = np.array([c * cb, c * (ca - cb * cg)/sg, v/(a * b * sg)])

        self._a_r = b * c * sa / v
        self._b_r = c * a * sb / v
        self._c_r = a * b * sg / v
        self._al_r = uacos((cb * cg - ca) / (sb * sg))
        self._be_r = uacos((cg * ca - cb) / (sg * sa))
        self._ga_r = uacos((ca * cb - cg) / (sa * sb))

        self._a_w = np.cross(self.b_v, self.c_v) / v
        self._b_w = np.cross(self.c_v, self.a_v) / v
        self._c_w = np.cross(self.a_v, self.b_v) / v
