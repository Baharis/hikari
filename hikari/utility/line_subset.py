def _min(*args):
    """Recursively returns minimum of arguments or their elements."""
    minimum = min(args)
    try:
        _ = iter(minimum)
    except TypeError:
        return minimum
    else:
        return _min(*minimum)


def _max(*args):
    """Recursively returns maximum of arguments or their elements."""
    maximum = max(args)
    try:
        _ = iter(maximum)
    except TypeError:
        return maximum
    else:
        return _max(*maximum)


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