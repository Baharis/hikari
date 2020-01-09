# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from collections import OrderedDict
from enum import Enum
import numpy as np
import pandas as pd
import pprint
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.umath import sin, cos
from uncertainties import unumpy as un
from kesshou.dataframes.cif import CifFrame
from kesshou.dataframes.hkl import HklCrystal


# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
"""This script compares two .cif files and calculates the difference between
their individual adp parameters as a correlation coefficient"""
# for now it works well only if 1 block is present

# Input details
input_cif1_name = 'PA_IAMref'
input_cif1_path = '/home/dtchon/x/HiPHAR/PA/IAM_models_simple/0kbar_full_experimental/0kbar_full_experimental.cif'
input_cif1_block = '0kbar_full_experimental'
input_cif2_name = 'ref'
input_cif2_path = '/home/dtchon/x/HiPHAR/PA/TAAM-shade_reference/plate_2_10.CIF'
input_cif2_block = 'I'

# Output details
output_directory = '/home/dtchon/x/HiPHAR/RFpirazB/IAM_models/'


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

def adp_det(u):
    """Calculate a determinant of an 3x3 matrix along with its uncertainty"""
    (u11, u12, u13) = (u[0, 0], u[0, 1], u[0, 2])
    (u21, u22, u23) = (u[1, 0], u[1, 1], u[1, 2])
    (u31, u32, u33) = (u[2, 0], u[2, 1], u[2, 2])
    det = u11*u22*u33 + u12*u23*u31 + u13*u21*u32 \
        - u31*u22*u13 - u32*u23*u11 - u33*u21*u12
    return det


def adp_cc(u, v):
    """Calculate cc of two ADP vectors based on doi:10.1107/S0907444999011853.
    Requires two 3x3 U matrices."""
    u = un.matrix(u)
    v = un.matrix(v)
    try:        #traktuj ujemne oraz zerowe jako całkowicie błędne
        if min(adp_det(u), adp_det(v)) <= 0:
            return ufloat(0., 0.)
        nominator = (adp_det(u.I) * adp_det(v.I)) ** (1 / 4)
        denominator = (adp_det(u.I + v.I) / 8) ** (1 / 2)
    except np.linalg.linalg.LinAlgError:
        return ufloat(0., 0.)
    return nominator / denominator


def adp_volume(u):
    """Take a 3x3 matrix and return a cube root of its determinant"""
    u = un.matrix(u)
    det = adp_det(u)
    try:
        sca = det ** (1 / 3)
    except ValueError:
        return ufloat(0., 0.)
    else:
        return sca


def adp_normalise(u):
    """Take a 3x3 matrix and normalise it to the determinant = 1"""
    u = un.matrix(u)
    try:
        reciprocal_volume = 1 / adp_volume(u)
    except ZeroDivisionError:
        reciprocal_volume = ufloat(0., 0.)
    for i in (0, 1, 2):
        for j in (0, 1, 2):
            u[i, j] *= reciprocal_volume
    return u


def adp_isotropise(u):
    """Take a 3x3 matrix and return a det*I matrix: det is input determinant"""
    u = un.matrix(u)
    return adp_volume(u) * un.matrix(np.identity(3))


class ReadingState(Enum):
    """An enumerator to describe cif reading state"""
    BEFORE = 1
    XYZ_KEYS = 2
    XYZ_VALUES = 3
    FILLER = 4
    ADP_KEYS = 5
    ADP_VALUES = 6
    AFTER = 7


def read_adps_from_cif(path):
    """Read adp and idp values from .cif file and output relevant dataframe"""
    lines = [line.strip().split() for line in open(path, 'r')]
    state = ReadingState.BEFORE
    xyz_labels = []
    u_iso = dict()
    data = OrderedDict([('_atom_site_aniso_label', list())])
    for line in lines:
        if state is ReadingState.BEFORE:
            try:
                if line[0] in ('_atom_site_label', '_atom_site_type_symbol'):
                    state = ReadingState.XYZ_KEYS
            except IndexError:
                pass

        if state is ReadingState.XYZ_KEYS:
            if line[0][0] == '_':
                xyz_labels.append(line[0])
                continue
            else:
                state = ReadingState.XYZ_VALUES

        if state is ReadingState.XYZ_VALUES:
            try:
                u_iso[line[xyz_labels.index('_atom_site_label')]] = \
                    line[xyz_labels.index('_atom_site_U_iso_or_equiv')]
            except (IndexError, ValueError):
                state = ReadingState.FILLER
                pass

        # if reading some unimportant garbage before adps
        if state is ReadingState.FILLER:
            try:
                if line[0] == '_atom_site_aniso_label':
                    state = ReadingState.ADP_KEYS
            except IndexError:
                pass
            continue

        # if reading _atom_site_aniso keys
        if state is ReadingState.ADP_KEYS:
            try:
                _ = ufloat_fromstr(line[1])
            except IndexError:
                data[line[0]] = []
            else:
                state = ReadingState.ADP_VALUES

        # if reading adp values
        if state is ReadingState.ADP_VALUES:
            try:
                _ = ufloat_fromstr(line[6])
            except IndexError:
                state = ReadingState.AFTER
            else:
                for index, key in enumerate(data.keys()):
                    if index not in (1, 2, 3, 4, 5, 6):
                        data[key].append(line[index])
                    else:
                        data[key].append(ufloat_fromstr(line[index]))

        # if reading some unimportant garbage after adps
        if state is ReadingState.AFTER:
            pass

    # fill the rest of missing adps with u_iso-scaled balls
    missing_data = list(set(u_iso.keys()) - set(data['_atom_site_aniso_label']))
    for missing in missing_data:
        zero = ufloat(0., 0.)
        u = ufloat_fromstr(u_iso[missing])
        if u.n < 0:
            scale = pow(u.n / abs(u.n) * u, 1 / 3)
        else:
            scale = pow(u, 1 / 3)
        data['_atom_site_aniso_label'].append(missing)
        data['_atom_site_aniso_U_11'].append(scale)
        data['_atom_site_aniso_U_22'].append(scale)
        data['_atom_site_aniso_U_33'].append(scale)
        data['_atom_site_aniso_U_12'].append(zero)
        data['_atom_site_aniso_U_23'].append(zero)
        data['_atom_site_aniso_U_13'].append(zero)

    return pd.DataFrame.from_dict(data=data, orient='columns')


# Read ADP information from the files
adps1 = read_adps_from_cif(input_cif1_path)
adps2 = read_adps_from_cif(input_cif2_path)


# Rename keys and sort the Pandas DataFrames'
def rename_and_sort_pandas_adp(adp):
    """Rename the pandas dataframe's columns to make the work easier later on"""
    renaming_dictionary = {
        '_atom_site_aniso_label': 'label',
        '_atom_site_aniso_U_11': 'xx',
        '_atom_site_aniso_U_22': 'yy',
        '_atom_site_aniso_U_33': 'zz',
        '_atom_site_aniso_U_12': 'xy',
        '_atom_site_aniso_U_13': 'xz',
        '_atom_site_aniso_U_23': 'yz',
    }
    adp.rename(index=str, columns=renaming_dictionary, inplace=True)
    adp.set_index('label', inplace=True)
    adp.sort_index(axis=0, inplace=True)
    return adp


adps1 = rename_and_sort_pandas_adp(adps1)
adps2 = rename_and_sort_pandas_adp(adps2)


def read_cell_from_cif(path, block):
    """use CifFrame to get cell parameters from the cif file"""

    # TODO exact angle precision is lowered (program assumes +/- 1, not exact)

    cif = CifFrame(file_path=path, file_data_block=block)
    a = ufloat_fromstr(cif.data['_cell_length_a'])
    b = ufloat_fromstr(cif.data['_cell_length_b'])
    c = ufloat_fromstr(cif.data['_cell_length_c'])
    al = ufloat_fromstr(cif.data['_cell_angle_alpha']) * np.pi / 180
    be = ufloat_fromstr(cif.data['_cell_angle_beta']) * np.pi / 180
    ga = ufloat_fromstr(cif.data['_cell_angle_gamma']) * np.pi / 180
    v = a * b * c * (1 - cos(al) ** 2 - cos(be) ** 2 - cos(ga) ** 2 +
                     2 * cos(al) * cos(be) * cos(ga)) ** 0.5
    return {'a': a, 'b': b, 'c': c, 'al': al, 'be': be, 'ga': ga, 'V': v}


# Transform fractional coordinates to cartesian coordinates:
def adp_fractional_to_cartesian_matrix(path_to_cif, block):
    """Transf. ADPs to cartesian based on doi:10.1107/S1600576715018075, A1"""
    cif = read_cell_from_cif(path_to_cif, block)
    a, b, c, v = cif['a'], cif['b'], cif['c'], cif['V']
    al, be, ga = cif['al'], cif['be'], cif['ga']
    m_11 = a
    m_12 = b * cos(ga)
    m_13 = c * cos(be)
    m_21 = 0
    m_22 = b * sin(ga)
    m_23 = c * (cos(al) - cos(be) * cos(ga)) / sin(ga)
    m_31 = 0
    m_32 = 0
    m_33 = v / (a * b * sin(ga))
    m = un.matrix([[m_11, m_12, m_13], [m_21, m_22, m_23], [m_31, m_32, m_33]])
    return m

# Create a lists of cc between the DataFrames:
def data_from_adps(adps1, adps2):
    """Calculate cc (correlation coeffitients = covariance) between adps"""

    # prepare objects and get transformation matrices
    data = {'label': [], 'cc': [], 'cc_size': [], 'cc_shape': [], 'scale': []}
    m1 = adp_fractional_to_cartesian_matrix(input_cif1_path, input_cif1_block)
    m2 = adp_fractional_to_cartesian_matrix(input_cif2_path, input_cif2_block)

    # prepare matrices and reorient them to cartesian system
    for (_, atom1), (_, atom2) in zip(adps1.iterrows(), adps2.iterrows()):
        assert atom1.name == atom2.name, 'Labels do not match'
        u1 = [[atom1['xx'], atom1['xy'], atom1['xz']],
              [atom1['xy'], atom1['yy'], atom1['yz']],
              [atom1['xz'], atom1['yz'], atom1['zz']]]
        u2 = [[atom2['xx'], atom2['xy'], atom2['xz']],
              [atom2['xy'], atom2['yy'], atom2['yz']],
              [atom2['xz'], atom2['yz'], atom2['zz']]]
        u1, u2 = un.matrix(u1), un.matrix(u2)
        u1 = u1 * m1
        u2 = u2 * m2
        # TODO - uncomment the following lines to get correct values
        # TODO it looks like coordinates change only increases uncertaint. (duh)

        # TODO - because the reorientation matrix is not normalised,
        # TODO comparing isotropic and unisotropic atom sizes GIVES NONSENSE

        # get the statistics
        data['label'].append(atom1.name)
        data['cc'].append(adp_cc(u=u1, v=u2))
        data['cc_size'].append(
            adp_cc(u=adp_isotropise(u1), v=adp_isotropise(u2)))
        data['cc_shape'].append(
            adp_cc(u=adp_normalise(u1), v=adp_normalise(u2)))
        #print(atom1.name, adp_volume(u1), adp_volume(u2))
        data['scale'].append(adp_volume(u1) / adp_volume(u2))

    data = pd.DataFrame.from_dict(data=data, orient='columns')
    data.set_index('label', inplace=True)
    return data


data = data_from_adps(adps1=adps1, adps2=adps2)


def data_divide(dataframe):
    """Return three dataframes: for H, non-H and all atoms"""
    hydrogens_indices, nonhydrogens_indices = [], []
    for label in dataframe.index:
        if label[0] == 'H' and not str.isalpha(label[1]):
            hydrogens_indices.append(label)
        else:
            nonhydrogens_indices.append(label)
    data_hydrogens = dataframe.drop(nonhydrogens_indices, axis=0)
    data_nonhydrogens = dataframe.drop(hydrogens_indices, axis=0)
    data_atoms = dataframe
    return data_hydrogens, data_nonhydrogens, data_atoms


def mean(iterable):
    """Arithmetic mean of elements present in iterable"""
    return (sum(iterable)) / max(len(iterable), 1)


def geomean(iterable):
    """Geometric mean of elements present in iterable"""
    value = 1.0
    for i in range(len(iterable)):
        value *= iterable[i]
    return pow(value, 1 / len(iterable))


def similarity(covariance):
    """Similarity index as defined by Whitten & Spackmann,
    for further details see F. Sanjuan et al., IURcJ, 3 (2016) 61-70"""
    return 100.0 * (1.0 - covariance)


def data_statistics(dataframe):
    """Calculate the "total" values of cc for investigated categories"""
    data_hydrogens, data_nonhydrogens, data_atoms = data_divide(dataframe)
    hydrogens = dict()
    hydrogens['cc'] = similarity(mean(data_hydrogens['cc']))
    hydrogens['cc_size'] = similarity(mean(data_hydrogens['cc_size']))
    hydrogens['cc_shape'] = similarity(mean(data_hydrogens['cc_shape']))
    hydrogens['scale'] = mean(data_hydrogens['scale'])
    dataframe = dataframe.append(
        pd.Series(data=hydrogens, name='S_hydrogens'))

    nonhydrogens = dict()
    nonhydrogens['cc'] = similarity(mean(data_nonhydrogens['cc']))
    nonhydrogens['cc_size'] = similarity(mean(data_nonhydrogens['cc_size']))
    nonhydrogens['cc_shape'] = similarity(mean(data_nonhydrogens['cc_shape']))
    nonhydrogens['scale'] = mean(data_nonhydrogens['scale'])
    dataframe = dataframe.append(
        pd.Series(data=nonhydrogens, name='S_nonhydrogens'))

    atoms = dict()
    atoms['cc'] = similarity(mean(data_atoms['cc']))
    atoms['cc_size'] = similarity(mean(data_atoms['cc_size']))
    atoms['cc_shape'] = similarity(mean(data_atoms['cc_shape']))
    atoms['scale'] = mean(data_atoms['scale'])
    dataframe = dataframe.append(
        pd.Series(data=atoms, name='S_atoms'))

    return dataframe


def data_to_csv(dataframe):
    """output .csv file with the ufloats divided into nominal and std_devs"""
    keys_n, keys_s = [], []
    for key in dataframe.keys():
        keys_n.append(key+'_n')
        keys_s.append(key+'_s')
    new_keys = ['label'] + keys_n + keys_s

    atoms = {'label': list()}
    old_keys = list(set(dataframe.keys()) - {'label'})
    for key in old_keys:
        atoms[key + '_n'] = list()
        atoms[key + '_s'] = list()
    for index, row in dataframe.iterrows():
        atoms['label'].append(row.name)
        for key in old_keys:
            atoms[key + '_n'].append(row[key].n)
            atoms[key + '_s'].append(row[key].s)
    atoms = pd.DataFrame.from_dict(data=atoms, orient='columns')
    atoms = atoms[new_keys]
    atoms.set_index('label', inplace=True)
    path = output_directory + input_cif1_name + '_vs_' + \
           input_cif2_name + '.csv'
    atoms.to_csv(path)


data = data_statistics(data)
data_to_csv(data)
print(data['cc_shape'])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
