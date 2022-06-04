"""
This file contains tools to work with tuples and lists used in the package.
"""

from numpy import power, linspace, sign, abs


def find_best(strings, criteria):
    """
    Parse a list of `strings` and return the "best" element based on `criteria`.

    :param strings: List of string where best will be found based on `criteria`.
    :type strings: List[str]
    :param criteria: '>'-separated substrings sought in descending order.
        '+' is a logical 'and', '=' substitutes: `A=B` returns B if A is found.
    :type criteria: str
    :return: Best string based on given `criteria` or `` if no best found.
    :rtype: str
    """
    criterion, _, further_criteria = criteria.partition('>')
    wanted = criterion.partition('=')[0].split('+')
    if all(w in strings or w is '' for w in wanted):
        return criterion.rpartition('=')[2]
    else:
        return find_best(strings, further_criteria)


def cubespace(start, stop=False, num=10, include_start=True):
    """
    Return sequence of *num* floats between *start* and *stop*.
    Analogously to numpy's linspace, values in returned list are chosen
    so that their cubes (hence name) are spread evenly in equal distance.

    If the parameter *stop* is not given, the value of *start* is used as
    the upper limit instead. In such case the lower limit is set to 0.

    The values of lower limit, *start*, and upper limit, *stop*,
    are included in the list. The *start* value can be excluded
    by setting the *include_start* keyword to False.

    :example:

    >>> cubespace(10, num=3)
    array([ 0.        ,  7.93700526,  10.       ])
    >>> cubespace(0, -10, num=3)
    array([ 0.        , -7.93700526, -10.       ])
    >>> cubespace(0, 10, num=3, include_start=False)
    array([ 6.93361274,  8.73580465, 10.        ])

    :param start: The starting value of a sequence.
    :type start: float
    :param stop: The ending value of a sequence. If False (default), *start*
        is used as a the ending value, while the starting value is set to 0.
    :type stop: float
    :param num: Number of samples to generate. Default is 10.
    :type num: int
    :param include_start: If False, the value of *start* is not included in
        returned list. Nonetheless, it is still considered as a starting point.
    :type include_start: bool
    :return: An array with *num* spaced samples in the *start*-*stop* interval.
    :rtype: numpy.ndarray
    """
    (start, stop) = (0.0, start) if stop is False else (start, stop)
    if include_start is False:
        return cubespace(start, stop, num=num+1, include_start=True)[1:]
    cubed_start = pow(start, 3)
    cubed_stop = pow(stop, 3)
    cubed_space = linspace(cubed_start, cubed_stop, num)
    return sign(cubed_space) * power(abs(cubed_space), 1/3)


def rescale_list_to_range(original, limits):
    """
    Linearly rescale values in original list to limits (minimum and maximum).

    :example:

    >>> rescale_list_to_range([1, 2, 3], (0, 10))
    [0.0, 5.0, 10.0]
    >>> rescale_list_to_range([1, 2, 3], (-10, 0))
    [-10.0, -5.0, 0.0]
    >>> rescale_list_to_range([1, 2, 3], (0j, 10j))
    [0j, 5j, 10j]

    :param original: Original list or list-like to be rescaled.
    :type original: list
    :param limits: Tuple of two floats, min and max, to constrain the new list
    :type limits: tuple
    :return: Original list rescaled to fit between min and max
    :rtype: list
    """
    new_min, new_max = limits[0:2]
    old_min, old_max = min(original), max(original)
    return (new_max + new_min) / 2 * original / old_min if old_min == old_max \
        else [new_max * (v - old_min) / (old_max - old_min) +
              new_min * (old_max - v) / (old_max - old_min) for v in original]


def rescale_list_to_other(source, target):
    """
    Linearly rescale *source* list of numeral values to
    elements of iterable scale in *target*.
    The numeric values in the first list are
    rescaled to the length of *other* using :func:`rescale_list_to_range`,
    changed to integers and used as pointers in *other* to retrieve final value.

    :example:

    >>> rescale_list_to_other([1, 2, 3], [-7.7, -6.6, -5.5, -4.4, -3.3])
    [-7.7, -5.5, -3.3]
    >>> rescale_list_to_other([-7.7, -6.6, -5.5, -4.4, -3.3], [1, 2, 3])
    [1, 1, 2, 3, 3]
    >>> rescale_list_to_other([1, 2, 3], 'holy grail')
    ['h', ' ', 'l']

    :param source: Iterable of numerals to be rescaled to *other*.
    :type source: Iterable
    :param target: Iterable from which the values in new list will be selected.
    :type target: Iterable
    :return: List with elements from *other* assigned using values in original.
    :rtype: list
    """
    return [target[int(v)] for v
            in rescale_list_to_range(source, (0, len(target)-1))]
