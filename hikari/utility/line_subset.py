def xmin(obj):
    try:
        return min(obj)
    except TypeError:
        return obj


def xmax(obj):
    try:
        return max(obj)
    except TypeError:
        return obj


class LineSubset:

    # SEGMENT METHODS
    class LineSegment:
        def __init__(self, left=None, right=None):
            if left >= right:
                raise ValueError('Left ({}) >= Right ({})'.format(left, right))
            self.left = left
            self.right = right

        def __eq__(self, other):
            return self.left == other.left and self.right == other.right

        def __ne__(self, other):
            return not(self == other)

        def __str__(self):
            return '[{}, {}]'.format(self.left, self.right)

        def __iter__(self):
            yield self.left
            yield self.right

        def __contains__(self, item):
            return all((xmin(self) <= xmin(item), xmax(self) >= xmax(item)))

    # CREATION METHODS
    def __init__(self, left=None, right=None):
        self.segments = []
        if left is None and right is None:
            pass
        elif right is None:
            for left1, left2 in left:
                self.segments.append(self.LineSegment(left1, left2))
        else:
            self.segments.append(self.LineSegment(left, right))
        self.segments = self._merge(self.segments)

    def _merge(self, segments):
        limits = sorted(list(zip(list(self), ['l', 'r'] * len(self.segments))))
        segments = []
        while limits:
            left = limits.pop(0)[0]
            depth = 1
            while limits:
                lim = limits.pop(0)
                if lim[1] == 'l':
                    depth += 1
                elif lim[1] == 'r' and depth > 1:
                    depth -= 1
                else:
                    segments.append(self.LineSegment(left, lim[0]))
                    break
        return segments

    # COMPARISON METHODS
    def __eq__(self, other):
        return self.segments == other.segments

    def __ne__(self, other):
        return not(self == other)

    def __lt__(self, other):
        return all(any(s in o for o in other.segments) for s in self.segments) \
               and not(self == other)

    def __gt__(self, other):
        return other < self

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return other < self or self == other

    # UNARY OPERATIONS
    def __pos__(self):
        return self

    def __neg__(self):
        limits = [-float('inf')] + list(self) + [float('inf')]
        lefts = limits[0::2]
        rights = limits[1::2]
        segments = [(lef, rig) for lef, rig in zip(lefts, rights) if lef < rig]
        return self.__class__(segments)

    # ARITHMETIC METHODS
    def __add__(self, other):
        segments = self.segments + other.segments
        return self.__class__(segments)

    def __sub__(self, other):
        return self * -other

    def __mul__(self, other):
        return -(-self + -other)

    def __truediv__(self, other):
        return (self - other) + (other - self)

    # REPRESENTATION METHODS
    def __str__(self):
        return ' + '.join([str(s) for s in self.segments])

    def __repr__(self):
        zipped = zip(list(self)[0::2], list(self)[1::2])
        pairs = ['{}, {}'.format(l, r) for l, r in zipped]
        return 'LineSubset({})'.format('; '.join(pairs))

    # CONTAINER METHODS
    def __iter__(self):
        return iter([lim for pair in self.segments for lim in pair])

    def __contains__(self, item):
        return any(item in segment for segment in self.segments)
