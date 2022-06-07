import numpy as np


class Interval:
    """
    A class handling [A, B] line segments and numpy operations on them.
    The following magic methods have been implemented for `Interval`:

    - creation using `Interval(A, B)` syntax (if `A == B`, represents a point)
    - comparison operations with numbers (True if holds for the whole Interval)
    - arithmetic operations, which act element-wise for left and right edge
    - access methods: `A == Interval(A, B)[0] == Interval(A, B).left`
    - container methods: 5 and 6.2 are `in Interval(5, 7)`, but 7. might not be
    """

    # CREATION METHODS
    def __init__(self, left, right):
        self.left = min(left, right)
        self.right = max(left, right)

    # COMPARISON METHODS
    def __eq__(self, other):
        return self.left == other.left and self.right == other.right

    def __lt__(self, other):
        return self.right < other

    def __gt__(self, other):
        return self.left > other

    def __le__(self, other):
        return self.right <= other

    def __ge__(self, other):
        return self.left >= other

    # UNARY OPERATIONS
    def __pos__(self):
        return self

    def __neg__(self):
        return Interval(-self.right, -self.left)

    # ARITHMETIC METHODS
    def __add__(self, other):
        return Interval(self.left + other, self.right + other)

    def __sub__(self, other):
        return Interval(self.left - other, self.right - other)

    def __mul__(self, other):
        return Interval(self.left * other, self.right * other)

    def __truediv__(self, other):
        return Interval(self.left / other, self.right / other)

    # REPRESENTATION METHODS
    def __str__(self):
        return f'[{self.left}, {self.right}]'

    def __repr__(self):
        return f'Interval({self.left}, {self.right})'

    # CONTAINER METHODS
    def __iter__(self):
        yield self.left
        yield self.right

    def __getitem__(self, key):
        return [self.left, self.right][key]

    def __setitem__(self, key, value):
        new_limits = [self.left, self.right]
        new_limits[key] = value
        self.left, self.right = new_limits

    def __contains__(self, item):
        return self.left <= self._min(item) and self._max(item) <= self.right

    def __len__(self):
        return 2

    @staticmethod
    def _min(item):
        """Recursively returns minimum of args if possible, otherwise args"""
        try:
            return min(item)
        except TypeError:
            return item

    @staticmethod
    def _max(item):
        """Recursively returns maximum of args if possible, otherwise args"""
        try:
            return max(item)
        except TypeError:
            return item

    def arange(self, step=1):
        """
        Return a 1D-list of values from left to right every step
        :param step: spacing between adjacent values, default 1.
        :type step: int or float
        :return: array of values from left to right (including right) every step
        :rtype: np.array
        """
        epsilon = step * 1e-7
        size = 1 + int((self.right - self.left + epsilon) / step)
        return np.array([self.left + step * i for i in range(size)])

    def comb_with(self, *others, step=1):
        """
        Return combinations of self.arange(step) with every other.arange(step)
        :param others: interval or iterable of intervals to comb self with
        :type others: Interval or tuple or list
        :param step: spacing between adjacent values, default 1.
        :type step: int or float
        :return: array of all combinations found in numpy.meshgrid every step
        :rtype: np.array
        """
        arrays = [self.arange(step)] + [other.arange(step) for other in others]
        return np.array(np.meshgrid(*arrays)).reshape(len(arrays), -1)

    def mesh_with(self, *others, step=1):
        """
        Return a numpy.mesh of self.arange(step) with every other.arange(step)
        :param others: interval or iterable of intervals to mesh self with
        :type others: Interval or tuple or list
        :param step: spacing between adjacent values, default 1.
        :type step: int or float
        :return: array of values meshed by numpy.meshgrid every step
        :rtype: np.array
        """
        arrays = [self.arange(step)] + [other.arange(step) for other in others]
        return np.meshgrid(*arrays)
