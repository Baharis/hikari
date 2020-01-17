import numpy as np
import numpy.linalg as lin
from .group import Group


class SpaceGroup(Group):
    """Basic Space Group class info holder"""

    unique_point = np.array([0.0 + np.pi / 100,
                             0.1 + np.pi / 100,
                             0.2 + np.pi / 100,
                             1])

    def __init__(self, generators):
        super().__init__(generators=generators)


# THIS CLASS IS NOT USED AT THE MOMENT; TO BE WORKED WITH LATER
