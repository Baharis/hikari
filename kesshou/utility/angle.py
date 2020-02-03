import numpy as np


def angle(value):
    """A simple function which interprets value in radian or degrees
    and returns value in radian. Interprets -2 < values < 2 as radians"""
    return value if -2 < value < 2 else np.radians(value)
