"""
This file contains class definition and necessary tools for space groups.
"""

import numpy as np
from .group import Group
from enum import Enum


class CrystalSystem(Enum):
    triclinic = 0
    monoclinic = 1
    orthorhombic = 2
    trigonal = 3
    tetragonal = 4
    cubic = 5
    hexagonal = 6


class SpaceGroup(Group):
    """
    Class containing information about three dimensional point groups. It
    inherits a lots of properties after :class:`kesshou.symmetry.group.Group`.
    At this moment this part of the code is not used and, due to high
    complexity, it is planned to introduce space groups in later versions.
    """

    unique_point = np.array([0.0 + np.pi / 100,
                             0.1 + np.pi / 100,
                             0.2 + np.pi / 100,
                             1])
    """Unique point in the space used for the sake of drawing etc."""

    def __init__(self, generators):
        super().__init__(generators=generators)
        operations = list()

    @property
    def system(self):
        folds = [op.fold for op in self.operations]
        orients = [op.orientation for op in self.operations]
        if 6 in folds:
            return CrystalSystem.hexagonal
        elif 3 in folds:
            orients_of_3 = len({o for f, o in zip(folds, orients) if f == 3})
            return CrystalSystem.cubic if orients_of_3 > 1 \
                else CrystalSystem.trigonal
        elif 4 in folds:
            return CrystalSystem.tetragonal
        elif 2 in folds:
            orients_of_2 = len({o for f, o in zip(folds, orients) if f == 2})
            return CrystalSystem.orthorhombic if orients_of_2 > 1 \
                else CrystalSystem.monoclinic
        else:
            return CrystalSystem.triclinic





# THIS CLASS IS NOT USED AT THE MOMENT; TO BE WORKED WITH LATER
