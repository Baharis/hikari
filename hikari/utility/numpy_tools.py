"""
This file contains basic tools to work with numpy objects used in the package.
"""

import re
import numpy as np


def str2array(s: str) -> np.ndarray:
    """
    Create a numpy int8 array using its compact representation saved as string.
    The input string might contain only digits and special characters.
    Every individual digit is treated as a new element, while the special
    characters induce special behaviour and include:

    - slash `/` - if present, divides main list into multiple sub-lists;
    - hyphen `-` - the next digit will be read as negative

    Please mind that the compact form used as this function's input is unfit
    to store floating point numbers or integers with absolute value above 9.

    :param s: input string containing only digits and special characters
    :return: a two-dimensional numpy array with integers
    """
    elements = [re.findall(r"-?\d", column) for column in s.split('/')]
    if len({len(element_column) for element_column in elements}) is not 1:
        raise ValueError('All columns specified with "/" must have same length')
    return np.array(elements, dtype=np.int8)
