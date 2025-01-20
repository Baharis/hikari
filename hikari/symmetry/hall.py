import re
from typing import List, Match, Union

import numpy as np

from hikari.symmetry.operations import BoundedOperation
from hikari.utility import dict_union


class HallSymbol:
    """Parse, interpret, and convert Hall symbol to `hikari.symmetry.Group`"""

    class HallSymbolException(Exception):
        """Exception raised when group can't be generated from the symbol"""

    REGEX = re.compile(
        r"""(-?)([pabcirstf])_(-?)(\d)([*xyz]?)([12345abcnuvwd]*)"""
        r"""(?:_(-?)(\d)(['"xyz]?)([12345abcnuvwd]*))?"""
        r"""(?:_(-?)(\d)\*?([12345abcnuvwd]*))?"""
        r"""(?:_-?(\d)([12345abcnuvwd]*))?(?:_\((\d)_(\d)_(\d)\))?""")
    r"""
    This is a regex which matches every possible lowercase Hall symbol for
    classical 3D space groups. It matches up to 17 groups in total as follows:

    1. Centrosymmetry: `[-]?`
    2. Lattice symbol: `[pabcirstf]` (required)
    3. Generator 1 inversion component: `[-]?`
    4. Generator 1 rotation fold: `\d` (required)
    5. Generator 1 rotation direction: `[*xyz]?`
    6. Generator 1 translation components: `[12345abcnuvwd]*`
    7. Generator 2 inversion component: `[-]?`
    8. Generator 2 rotation fold: `\d`
    9. Generator 2 rotation direction: `['"xyz]?`
    10. Generator 2 translation components: `[12345abcnuvwd]*`
    11. Generator 3 inversion component: `[-]?`
    12. Generator 3 rotation fold: `\d`
    13. Generator 3 translation components: `[12345abcnuvwd]*`
    14. Generator 4 rotation fold: `\d`
    15. Generator 4 translation components: `[12345abcnuvwd]*`
    16. Origin shift constituent 1 expressed as count of 1/12 shifts: `\d`
    17. Origin shift constituent 2 expressed as count of 1/12 shifts: `\d`
    18. Origin shift constituent 3 expressed as count of 1/12 shifts: `\d`
    """

    DIRECTION_SYMBOLS = ("'", '"', '*')
    TRANSLATION_SYMBOLS = '12345abcnuvwd'
    LATTICE_GENERATORS = {
        'p': (BoundedOperation(np.eye(3)),),
        'a': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (0, 1 / 2, 1 / 2))),
        'b': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (1 / 2, 0, 1 / 2))),
        'c': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (1 / 2, 1 / 2, 0))),
        'i': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (1 / 2, 1 / 2, 1 / 2))),
        'r': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (2 / 3, 1 / 3, 1 / 3)),
              BoundedOperation(np.eye(3), (1 / 3, 2 / 3, 2 / 3))),
        's': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (1 / 3, 1 / 3, 2 / 3)),
              BoundedOperation(np.eye(3), (2 / 3, 2 / 3, 1 / 3))),
        't': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (1 / 3, 2 / 3, 1 / 3)),
              BoundedOperation(np.eye(3), (2 / 3, 1 / 3, 2 / 3))),
        'f': (BoundedOperation(np.eye(3)),
              BoundedOperation(np.eye(3), (0, 1 / 2, 1 / 2)),
              BoundedOperation(np.eye(3), (1 / 2, 0, 1 / 2)),
              BoundedOperation(np.eye(3), (1 / 2, 1 / 2, 0))),
    }
    UNIVERSAL_MATRICES = {
        '1': np.eye(3),
    }
    PRINCIPAL_ROTATIONS = {
        'x': {'1': np.eye(3),
              '2': np.diag([1, -1, -1]),
              '3': np.array([(1, 0, 0), (0, 0, -1), (0, 1, -1)]),
              '4': np.array([(1, 0, 0), (0, 0, -1), (0, 1, 0)]),
              '6': np.array([(1, 0, 0), (0, 1, -1), (0, 1, 0)]),},
        'y': {'1': np.eye(3),
              '2': np.diag([-1, 1, -1]),
              '3': np.array([(-1, 0, 1), (0, 1, 0), (-1, 0, 0)]),
              '4': np.array([(0, 0, 1), (0, 1, 0), (-1, 0, 0)]),
              '6': np.array([(0, 0, 1), (0, 1, 0), (-1, 0, 1)]),},
        'z': {'1': np.eye(3),
              '2': np.diag([-1, -1, 1]),
              '3': np.array([(0, -1, 0), (1, -1, 0), (0, 0, 1)]),
              '4': np.array([(0, -1, 0), (1, 0, 0), (0, 0, 1)]),
              '6': np.array([(1, -1, 0), (1, 0, 0), (0, 0, 1)]),},
        '*': {'3': np.array([(0, 0, 1), (1, 0, 0), (0, 1, 0)]),}
    }
    FACE_X_DIAGONAL_ROTATIONS = {
        "1'": np.eye(3),
        '2"': np.array([(-1, 0, 0), (0, 0, 1), (0, 1, 0)]),
        "2'": np.array([(-1, 0, 0), (0, 0, -1), (0, -1, 0)]),
    }
    FACE_Y_DIAGONAL_ROTATIONS = {
        "1'": np.eye(3),
        '2"': np.array([(0, 0, 1), (0, -1, 0), (1, 0, 0)]),
        "2'": np.array([(0, 0, -1), (0, -1, 0), (-1, 0, 0)]),
    }
    FACE_Z_DIAGONAL_ROTATIONS = {
        "1'": np.eye(3),
        '2"': np.array([(0, 1, 0), (1, 0, 0), (0, 0, -1)]),
        "2'": np.array([(0, -1, 0), (-1, 0, 0), (0, 0, -1)]),
    }
    BODY_DIAGONAL_ROTATIONS = {
        '1': np.eye(3),
        '3': np.array([(0, 0, 1), (1, 0, 0), (0, 1, 0)]),
    }
    STATIC_TRANSLATIONS = {
        'a': np.array([1/2, 0, 0]),
        'b': np.array([0, 1/2, 0]),
        'c': np.array([0, 0, 1/2]),
        'n': np.array([1/2, 1/2, 1/2]),
        'u': np.array([1/4, 0, 0]),
        'v': np.array([0, 1/4, 0]),
        'w': np.array([0, 0, 1/4]),
        'd': np.array([1 / 4, 1 / 4, 1 / 4]),
    }
    DYNAMIC_TRANSLATIONS = {
        '3x': {'1': np.array([1 / 3, 0, 0]), '2': np.array([2 / 3, 0, 0])},
        '3y': {'1': np.array([0, 1 / 3, 0]), '2': np.array([0, 2 / 3, 0])},
        '3z': {'1': np.array([0, 0, 1 / 3]), '2': np.array([0, 0, 2 / 3])},
        '4x': {'1': np.array([1 / 4, 0, 0]), '3': np.array([3 / 4, 0, 0])},
        '4y': {'1': np.array([0, 1 / 4, 0]), '3': np.array([0, 3 / 4, 0])},
        '4z': {'1': np.array([0, 0, 1 / 4]), '3': np.array([0, 0, 3 / 4])},
        '6x': {'1': np.array([1 / 6, 0, 0]), '2': np.array([1 / 3, 0, 0]),
               '4': np.array([2 / 3, 0, 0]), '5': np.array([5 / 6, 0, 0])},
        '6y': {'1': np.array([0, 1 / 6, 0]), '2': np.array([0, 1 / 3, 0]),
               '4': np.array([0, 2 / 3, 0]), '5': np.array([0, 5 / 6, 0])},
        '6z': {'1': np.array([0, 0, 1 / 6]), '2': np.array([0, 0, 1 / 3]),
               '4': np.array([0, 0, 2 / 3]), '5': np.array([0, 0, 5 / 6])},
    }

    def __init__(self, hall_symbol: str) -> None:
        self.symbol = hall_symbol

    @property
    def elements(self) -> Union[None, Match]:
        """Return `self.symbol` elements indexed as in `self.HALL_REGEX` doc"""
        return self.REGEX.match(self.symbol)

    @property
    def symbol(self) -> str:
        return self._symbol

    @symbol.setter
    def symbol(self, symbol: str) -> None:
        self._symbol = symbol.lower().replace(' ', '_')

    @property
    def generators(self) -> List[BoundedOperation]:
        elements: re.Match = self.elements

        # lattice / centering generators
        generators: List[BoundedOperation] = list(self.LATTICE_GENERATORS[elements[2]])
        if elements[1] == '-':
            generators.append(BoundedOperation(-np.eye(3)))

        # generator 1
        gen1_inv, gen1_fold, gen1_dir, gen1_tl = elements.group(3, 4, 5, 6)
        tf1 = -np.eye(3) if gen1_inv else np.eye(3)
        gen1_dir = gen1_dir or 'z'
        tf1 = tf1 @ self.PRINCIPAL_ROTATIONS[gen1_dir][gen1_fold]
        tl1 = np.array([0., 0., 0.])
        gen1_translation_dicts = [self.STATIC_TRANSLATIONS]
        if gen1_fold in '346':
            screw_tls = self.DYNAMIC_TRANSLATIONS.get(gen1_fold + gen1_dir, {})
            gen1_translation_dicts.append(screw_tls)
        gen1_translations = dict_union(*gen1_translation_dicts)
        for tl in gen1_tl:
            tl1 += gen1_translations[tl]
        generators.append(BoundedOperation(tf1, tl1))

        # generator 2
        gen2_inv, gen2_fold, gen2_dir, gen2_tl = elements.group(7, 8, 9, 10)
        if gen2_fold:
            tf2 = -np.eye(3) if gen2_inv else np.eye(3)
            if gen1_fold in '24':
                gen2_dir = gen2_dir or "x"
                tf2 = tf2 @ self.PRINCIPAL_ROTATIONS[gen2_dir][gen2_fold]
            else:  # if gen1_fold in '36'
                gen2_dir = gen2_dir or "'"
                diagonal_rots = self.FACE_X_DIAGONAL_ROTATIONS if gen1_dir == 'x' \
                    else self.FACE_Y_DIAGONAL_ROTATIONS if gen1_dir == 'y' \
                    else self.FACE_Z_DIAGONAL_ROTATIONS
                tf2 = tf2 @ diagonal_rots[gen2_fold + gen2_dir]
            tl2 = np.array([0., 0., 0.])
            for tl in gen2_tl:
                tl2 += self.STATIC_TRANSLATIONS[tl]
            generators.append(BoundedOperation(tf2, tl2))

        # generator 3
        gen3_inv, gen3_fold, gen3_tl = elements.group(11, 12, 13)
        if gen3_fold:
            tf3 = -np.eye(3) if gen3_inv else np.eye(3)
            tf3 = tf3 @ self.BODY_DIAGONAL_ROTATIONS[gen3_fold]
            tl3 = np.array([0., 0., 0.])
            for tl in gen3_tl:
                tl3 += self.STATIC_TRANSLATIONS[tl]
            generators.append(BoundedOperation(tf3, tl3))

        # generator 4
        gen4_fold, gen4_tl = elements.group(14, 15)
        if gen4_fold:
            tl4 = np.array([0., 0., 0.])
            for tl in gen4_tl:
                tl4 += self.STATIC_TRANSLATIONS[tl]
            generators.append(BoundedOperation(-np.eye(3), tl4))

        if any(elements.group(16, 17, 18)):
            shift_ints = [int(i) for i in elements.group(16, 17, 18)]
            shift_vector = np.array(shift_ints, dtype=float) / 12
            shift_op_left = BoundedOperation(np.eye(3), shift_vector)
            shift_op_right = BoundedOperation(np.eye(3), -shift_vector)
            generators = [BoundedOperation.from_matrix(shift_op_left.matrix @ g.matrix
                                                @ shift_op_right.matrix) for g in generators]

        return generators
