"""
This file contains constants and functions containing operating system and path
information / methods used in the package.
"""

import os


def make_absolute_path(*path_elements):
    """
    Return an absolute path to a specified file.
    Can be supplied a list of elements to join.
    In case of receiving a single string, an absolute path to the location
    pointed by the string will be returned. In case of receiving multiple
    arguments, they will be joined using system-specific syntax beforehand.

    :Example:

    >>> make_absolute_path('')
    '/home/user/hikari'
    >>> make_absolute_path('..', 'flying_circus', 'loc.txt')
    '/home/user/flying_circus/loc.txt'
    >>> make_absolute_path('~')
    '/home/user/hikari/~'

    :param path_elements: path element (string or path-like) or multiple
        path elements to be joined together
    :type path_elements: str
    :return: absolute path to the given location
    :rtype: str
    """
    return os.path.abspath(os.path.join(*path_elements))


home_directory = os.path.expanduser('~')
