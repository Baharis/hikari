"""
This file contains all constants and functions containing purely chemical
information / methods used in the package.
"""

import re

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


def split_atom_label(label):
    """
    Split atomic label used in ref or cif formats into computer-readable parts.
    :param label: Short atomic label such as "C12A" or "Fe1".
    :type label: str
    :return: 3-element tuple with element, number, and suffix, as str.
    :rtype: tuple
    """
    label_regex = re.compile(r'(\d+|\s+)')
    label_fragments = label_regex.split(label)
    element = label_fragments[0]
    number = label_fragments[1]
    suffix = ''.join(label_fragments[2:])
    return element, number, suffix
