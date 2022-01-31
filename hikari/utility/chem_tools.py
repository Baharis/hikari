"""
This file contains all constants and functions containing purely chemical
information / methods used in the package.
"""
from numpy import exp as e
from hikari.resources import Xray_atomic_form_factors

chemical_elements = \
    ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
     'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
     'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
     'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og')
"""Tuple containing 2-letter symbols of the first 118 chemical elements."""


def atomic_form_factor(atom, r):
    f = Xray_atomic_form_factors.loc[atom]
    s2 = r**2 / 4
    return f['a1'] * e(-f['b1'] * s2) + f['a2'] * e(-f['b2'] * s2) + \
           f['a3'] * e(-f['b3'] * s2) + f['a4'] * e(-f['b4'] * s2) + f['c']
