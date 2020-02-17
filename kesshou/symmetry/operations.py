"""
This file contains the definition of 3-dimensional point symmetry operations.
"""

import numpy as np


symm_ops = dict()
"""
Dictionary containing 3x3 matrices of all 3-dimensional point-symmetry 
generators. Please mind that this file does not contain all symmetry operations,
but only those that are necessary to generate the point groups.

The names of symmetry operations follow the international crystallographic
naming convention wherever possible. For trigonal and hexagonal system case,
the names are additionally prefixed with the letter 'h'. 
"""
symm_ops['1'] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
symm_ops['-1'] = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
symm_ops['2_x'] = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
symm_ops['2_y'] = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
symm_ops['2_z'] = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
symm_ops['2_xy'] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
symm_ops['-3_z'] = np.array([[1, 1, 0], [-1, 0, 0], [0, 0, -1]])
symm_ops['3_xyz'] = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
symm_ops['-3_xyz'] = np.array([[0, 0, -1], [-1, 0, 0], [0, -1, 0]])
symm_ops['4_z'] = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
symm_ops['-4_z'] = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
symm_ops['m_x'] = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
symm_ops['m_y'] = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
symm_ops['m_z'] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
symm_ops['m_xy'] = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
symm_ops['h2_x'] = np.array([[1, -1, 0], [0, -1, 0], [0, 0, -1]])
symm_ops['h2_x2y'] = np.array([[-1, 1, 0], [0, 1, 0], [0, 0, -1]])
symm_ops['h3_z'] = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
symm_ops['h-3_z'] = np.array([[0, 1, 0], [-1, 1, 0], [0, 0, -1]])
symm_ops['h6_z'] = np.array([[1, -1, 0], [1, 0, 0], [0, 0, 1]])
symm_ops['h-6_z'] = np.array([[-1, 1, 0], [-1, 0, 0], [0, 0, -1]])
symm_ops['hm_x'] = np.array([[-1, 1, 0], [0, 1, 0], [0, 0, 1]])
symm_ops['hm_x2y'] = np.array([[1, -1, 0], [0, -1, 0], [0, 0, 1]])
symm_ops['hm_z'] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
