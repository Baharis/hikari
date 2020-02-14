import os


def make_absolute_path(*path_elements):
    return os.path.abspath(os.path.join(*path_elements))


home_directory = os.path.expanduser('~')
