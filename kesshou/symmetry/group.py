import numpy as np
import numpy.linalg as lin
from itertools import product as itertools_product


class Group:
    """Basic Group class info holder"""

    unique_point = np.array([1])

    def __init__(self, generators):
        self.generators = generators
        self.operations = list()
        self.equivalents = list()
        self.construct()

    @property
    def chiral_operations(self):
        return [op for op in self.operations if lin.det(op) > 0]

    def construct(self):
        self.generate_operations()
        self.generate_equivalents()

    def generate_operations(self):
        """generate whole point group based on its generators"""
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

    def generate_equivalents(self):
        """generate all equivalent points based on its operations"""
        eqs = list()
        for op in self.operations:
            eqs.append(np.dot(op, self.unique_point))
        self.equivalents = eqs

    @property
    def is_chiral(self):
        return all(lin.det(op) > 0 for op in self.operations)

    @property
    def is_polar(self):
        sum_of_equivalents = sum(self.equivalents)
        return lin.norm(sum_of_equivalents) > 0.1
