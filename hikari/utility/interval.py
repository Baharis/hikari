def _min(*args):
    """Recursively returns minimum of args if possible, otherwise args"""
    try:
        return min(args)
    except TypeError:
        return args


def _max(*args):
    """Recursively returns maximum of args if possible, otherwise args"""
    try:
        return max(args)
    except TypeError:
        return args


class Interval:
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

    def __contains__(self, item):
        return self.left <= _min(item) and _max(item) <= self.right

    def arange(self, step=1):
        """
        Return a 1D-list of values from left to right every step
        :param step: spacing between adjacent values, default 1.
        :type step: int or float
        :return: list of values from left to right (including right) every step
        :rtype: list
        """
        epsilon = step * 1e-7
        size = 1 + int((self.right - self.left + epsilon) / step)
        return [self.left + step * i for i in range(size)]









