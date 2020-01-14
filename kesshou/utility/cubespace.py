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

