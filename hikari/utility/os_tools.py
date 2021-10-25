"""
This file contains constants and functions containing operating system and path
information / methods used in the package.
"""

import os, pathlib


def make_abspath(*path_elements):
    """
    Return an absolute path to a specified file. Accepts one or more strings,
    which are joined according to system-specific syntax. Relative paths are
    intepreted as having root at current working directory. Links and symbols
    such as '~', '..' or '/' are expanded and interpreted. Function calls with
    nothing and '~' as arguments yield working and home directory, respectively.

    :Example:

    >>> make_abspath()
    PosixPath('/home/username/hikari/hikari/utility')
    >>> make_abspath('~')
    PosixPath('/home/username')
    >>> make_abspath('~', 'cocoa/is', 'not/../hot.txt')
    PosixPath('/home/username/cocoa/is/hot.txt')

    Running the script on Windows machine will yield an appropriate Windows
    Path. In particular, 'cocoa/is' will be intepreted as single directory.

    :param path_elements: path element(s) to join and intepret, defaults to cwd.
    :type path_elements: str
    :return: absolute path to the given location
    :rtype: pathlib.Path
    """
    return pathlib.Path().joinpath(*path_elements).expanduser().resolve()


home_directory = os.path.expanduser('~')
