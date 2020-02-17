"""
This file contains class definition and necessary tools for space groups.
"""

import numpy as np
import numpy.linalg as lin
from .group import Group


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


# THIS CLASS IS NOT USED AT THE MOMENT; TO BE WORKED WITH LATER
