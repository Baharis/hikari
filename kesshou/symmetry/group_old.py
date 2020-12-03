"""
This file contains class definition and necessary tools for constructing
and evaluating all symmetry groups.
"""
import numpy as np
import numpy.linalg as lin
from itertools import product as itertools_product


class Group:
    """
    Base class containing information about symmetry groups. It acts as a base
    object for :class:`kesshou.symmetry.pointgroup.PointGroup` and
    :class:`kesshou.symmetry.spacegroup.SpaceGroup`.

    In the future it is planned to expand this group using getter/setter
    operations to automatically refresh operations and symmetry equivalents
    whenever a list of generators is changed.
    """

    unique_point = np.array([1])
    """Unique point to be transformed by symmetry operations."""

    def __init__(self, generators):
        self.generators = generators
        """A list of all group generators."""
        self.operations = list()
        """A list of all group elements."""
        self.equivalents = list()
        """A list of all points equivalent by symmetry."""
        self.construct()

    @property
    def chiral_operations(self):
        """
        A list of symmetry operations which preserve structure chirality.

        :return: Symmetry operations whose matrix determinant is positive.
        :rtype: list
        """
        return [op for op in self.operations if lin.det(op) > 0]

    def construct(self):
        """
        Prepare a list of symmetry group operations and symmetry equivalent
        points using the list of generators.
        """
        self._generate_operations()
        self._generate_equivalents()

    def _generate_operations(self):
        """Generate a list of operations based on the generators."""
        operations = self.generators

        def _is_in(op, ops):
            return any([np.array_equal(op, op1) for op1 in ops])

        # then perform recursive group table multiplication
        def _find_new_product(current_ops):
            new_op = None
            for op1, op2 in itertools_product(current_ops, current_ops):
                op = np.dot(op1, op2)
                if not _is_in(op, current_ops):
                    new_op = op
                    break
            if new_op is None:
                return current_ops
            else:
                return _find_new_product([*current_ops, new_op])

        operations = _find_new_product(operations)
        self.operations = operations

    def _generate_equivalents(self):
        """Generate a list of equivalent points based on group operations."""
        eqs = list()
        for op in self.operations:
            eqs.append(np.dot(op, self.unique_point))
        self.equivalents = eqs

    @property
    def is_chiral(self):
        """
        Check whether all group operations preserve the chirality.

        :return: True if matrix determinant of all group operations is positive,
            False otherwise.
        :rtype: bool
        """
        return all(lin.det(op) > 0 for op in self.operations)

    @property
    def is_polar(self):
        """
        Check whether all group operations preserve the polarity.

        :return: True if a general polar property is preserved,
            e.g. if a sum of symmetry equivalent vectors is 0, False otherwise.
        :rtype: bool
        """
        sum_of_equivalents = sum(self.equivalents)
        return lin.norm(sum_of_equivalents) > 0.1

# TODO Think about integrating / using / following "xcore" package -
# TODO very useful, but must be rewritten from Python2/c
