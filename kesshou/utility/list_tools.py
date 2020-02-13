from numpy import power, linspace


def cubespace(start, stop=False, num=10, include_start=True):
    """Return num floats between start and stop (including start and stop)
    whose cubic roots are evenly spaced"""

    # change start to 0 if only one value is given
    if stop is False:
        stop = start
        start = 0.0

    # swap positions of start and stop if start > stop
    start = min(start, stop)
    stop = max(start, stop)

    # add one point to num if start is not to be included
    num = num+1 if not include_start else num

    # take the cube values of start and stop
    start3 = pow(start, 3)
    stop3 = pow(stop, 3)

    # obtain 'num' values from 'start3' to 'stop3'
    space3 = linspace(start3, stop3, num)

    # get the cubic root of all values
    space = power(space3, 1/3)

    # remove start if it should not be included
    space = space[1:] if not include_start else space

    # return obtained points
    return space


def rescale_list_to_range(original, limits):
    """
    Linearly rescale values in original list to limits (minimum and maximum).

    :param original: Original list or list-like to be rescaled.
    :type original: list
    :param limits: Tuple of two floats, min and max, to constrain the new list
    :type limits: tuple
    :return: Original list rescaled to fit between minimum and maximum
    :rtype: list
    """
    new_min, new_max = limits[0:2]
    old_min, old_max = min(original), max(original)
    return [new_max * (v - old_min) / (old_max - old_min) +
            new_min * (old_max - v) / (old_max - old_min) for v in original]


def rescale_list_to_other(original, other):
    """
    Linearly rescale original list to a template scale of elements in other.

    :param original: List to be rescaled to fixed scale.
    :type original: list
    :param other: List of ordered values to act as a scale according to which
        original list is rescaled based on its values
    :type other: list
    :return: List with elements from scale assigned using values in original.
    :rtype: list
    """
    minimum, maximum = min(original), 0.99999999 * max(original)
    if maximum <= minimum:
        return other[:1] * len(original)
    else:
        return [other[int(val)] for val in
                rescale_list_to_range(original, (0, 0.99999 * len(other)))]
