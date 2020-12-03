from kesshou.symmetry import Group, SymmOp
import numpy as np

# short notation to effectively re-use some symmetry elements in group creation
# r & ri denote rotation and rotoinvertion in monoclinic or orthogonal systems
# h & hi denote rotation and rotoinvertion in hexagonal system only
# CENTERING
A = SymmOp.from_code('x, y+1/2, z+1/2')
B = SymmOp.from_code('x+1/2, y, z+1/2')
C = SymmOp.from_code('x+1/2, y+1/2, z')
I = SymmOp.from_code('x+1/2, y+1/2, z+1/2')
# ROTATIONS
r2y = SymmOp.from_code('-x, y, -z')
r21y = SymmOp.from_code('-x, y+1/2, -z')
r2z = SymmOp.from_code('-x, -y, z')
r21z = SymmOp.from_code('-x, -y, z+1/2')
r4z = SymmOp.from_code('-y, x, z')
r41z = SymmOp.from_code('-y, x, z+1/4')
h3z = SymmOp.from_code('-y, x-y, z')
h31z = SymmOp.from_code('-y, x-y, z+1/3')
h6z = SymmOp.from_code('x-y, x, z')
h61z = SymmOp.from_code('x-y, x, z+1/6')
# REFLECTIONS
ay = SymmOp.from_code('x+1/2, -y, z')
bx = SymmOp.from_code('-x, y+1/2, z')
cx = SymmOp.from_code('-x, y, z+1/2')
cy = SymmOp.from_code('x, -y, z+1/2')
dx = SymmOp.from_code('-x, y+1/4, z+1/4')
mx = SymmOp.from_code('-x, y, z')
my = SymmOp.from_code('x, -y, z')
mz = SymmOp.from_code('x, y, -z')
ny = SymmOp.from_code('x+1/2, -y, z+1/2')
# INVERSION
i = SymmOp.from_code('-x,-y,-z')
# DIRECTIONS
x = np.array((1/4, 0, 0))
y = np.array((0, 1/4, 0))
z = np.array((0, 0, 1/4))

SG = {
    # TRICLINIC
    'P1': Group(SymmOp.from_code('x,y,z')),
    'P-1': Group(i),
    # MONOCLINIC
    'P2': Group(r2y),
    'P21': Group(r21y),
    'C2': Group(C, r2y),
    'Pm': Group(my),
    'Pc': Group(cy),
    'Cm': Group(C, my),
    'Cc': Group(C, cy),
    'P2/m': Group(r2y, my),
    'P21/m': Group(r21y, my),
    'C2/m': Group(C, r2y, my),
    'P2/c': Group(r2y.at(z), cy),
    'P21/c': Group(r21y.at(z), cy),
    'C2/c': Group(C, r2y.at(z), cy),
    # ORTHORHOMBIC
    'P222': Group(r2z, r2y),
    'P2221': Group(r21z, r2y.at(z)),
    'P21212': Group(r2z, r21y.at(x)),
    'P212121': Group(r21z.at(x), r21y.at(z)),
    'C2221': Group(C, r21z, r2y.at(z)),
    'C222': Group(C, r2z, r2y),
    'F222': Group(A, B, r2z, r2y),
    'I222': Group(I, r2z, r2y),
    'I212121': Group(I, r2z.at(y), r2y.at(x)),
    'Pmm2': Group(r2z, mx),
    'Pmc21': Group(r21z, mx),
    'Pcc2': Group(r2z, cy),
    'Pma2': Group(r2z, mx.at(x)),
    'Pca21': Group(r21z, cx.at(x)),
    'Pnc2': Group(r2z, cy.at(y)), #last
    'Pmn21': Group(r21z.at(x), mx),
    'Pba2': Group(r2z, bx.at(x)),
    'Pna21': Group(r21z, ay.at(y)),
    'Pnn2': Group(r2z, ny.at(y)),
    'Cmm2': Group(C, r2z, mx),
    'Cmc21': Group(C, r21z, mx),
    'Ccc2': Group(C, r2z, cx),
    'Amm2': Group(A, r2z, mx),
    'Aem2': Group(A, r2z, cy),
    'Ama2': Group(A, r2z, ay),
    'Aea2': Group(A, r2z, ay.at(y)),
    'Fmm2': Group(A, B, r2z, mx),
    'Fdd2': Group(A, B, r2z, dx.at(np.array([1/8, 0, 0]))), #
    'Imm2': Group(I, r2z, mx),
    'Iba2': Group(I, r2z, cx),
    'Ima2': Group(I, r2z, ay),
    # TETRAGONAL
    # TRIGONAL
    # HEXAGONAL
    # CUBIC
}
"""
Dictionary containing all known space groups written as :class:`Group`
along with alternative axis settings. The point groups in this dictionary
can be accessed using their short Hermann-Maugin notation, as presented below.
The origin is always as suggested by http://img.chem.ucl.ac.uk/sgp/large/sgp.htm

+-------+--------------+----------------+----------------+-----------------+
| No.   | CRYSTAL      | Hermann-Maugin | Hermann-Maugin | Can be accessed |
|       | SYSTEM       | short notation | full notation  | using           |
+-------+--------------+----------------+----------------+-----------------+
| 1     | triclinic    | P1             | P 1            | `SG['P1']`      |
+-------+              +----------------+----------------+-----------------+
| 2     |              | P-1            | P -1           | `SG['P-1']`     |
+-------+--------------+----------------+----------------+-----------------+
| 3     | monoclinic   | P2             | P 1 2 1        | `SG['P2']`      |
+-------+              +----------------+----------------+-----------------+
| 4     |              | P21            | P 1 21 1       | `SG['P21']`     |
+-------+              +----------------+----------------+-----------------+
| 5     |              | C2             | C 1 2 1        | `SG['C2']`      |
+-------+              +----------------+----------------+-----------------+
| 6     |              | Pm             | P 1 m 1        | `SG['Pm']`      |
+-------+              +----------------+----------------+-----------------+
| 7     |              | Pc             | P 1 c 1        | `SG['Pc']`      |
+-------+              +----------------+----------------+-----------------+
| 8     |              | Cm             | C 1 m 1        | `SG['Cm']`      |
+-------+              +----------------+----------------+-----------------+
| 9     |              | Cc             | C 1 c 1        | `SG['Cc']`      |
+-------+              +----------------+----------------+-----------------+
| 10    |              | P2/m           | P 1 2/m 1      | `SG['P2/m']`    |
+-------+              +----------------+----------------+-----------------+
| 11    |              | P21/m          | P 1 21/m 1     | `SG['P21/m']`   |
+-------+              +----------------+----------------+-----------------+
| 12    |              | C2/m           | C 1 2/m 1      | `SG['C2/m']`    |
+-------+              +----------------+----------------+-----------------+
| 13    |              | P2/c           | P 1 2/c 1      | `SG['P2/c']`    |
+-------+              +----------------+----------------+-----------------+
| 14    |              | P21/c          | P 1 21/c 1     | `SG['P21/c']`   |
+-------+              +----------------+----------------+-----------------+
| 15    |              | C2/c           | C 1 2/c 1      | `SG['C2/c']`    | 
+-------+--------------+----------------+----------------+-----------------+
| 16    | orthorhombic | P222           | P 2 2 2        | `SG['P222']`    |
+-------+              +----------------+----------------+-----------------+
| 17    |              | P2221          | P 2 2 21       | `SG['P2221']`   |
+-------+              +----------------+----------------+-----------------+
| 18    |              | P21212         | P 21 21 2      | `SG['P21212']`  |
+-------+              +----------------+----------------+-----------------+
| 19    |              | P212121        | P 21 21 21     | `SG['P212121']` |
+-------+              +----------------+----------------+-----------------+
| 20    |              | C2221          | C 2 2 21       | `SG['C2221']`   |
+-------+              +----------------+----------------+-----------------+
| 21    |              | C222           | C 2 2 2        | `SG['C222']`    |
+-------+              +----------------+----------------+-----------------+
| 22    |              | F222           | F 2 2 2        | `SG['F222']`    |
+-------+              +----------------+----------------+-----------------+
| 23    |              | I222           | I 2 2 2        | `SG['I222']`    |
+-------+              +----------------+----------------+-----------------+
| 24    |              | I212121        | I 21 21 21     | `SG['I212121']` |
+-------+              +----------------+----------------+-----------------+
| 25    |              | Pmm2           | P m m 2        | `SG['Pmm2']`    |
+-------+              +----------------+----------------+-----------------+
| 26    |              | Pmc21          | P m c 21       | `SG['Pmc21']`   |
+-------+              +----------------+----------------+-----------------+
| 27    |              | Pcc2           | P c c 2        | `SG['Pcc2']`    |
+-------+              +----------------+----------------+-----------------+
| 28    |              | Pma2           | P m a 2        | `SG['Pma2']`    |
+-------+              +----------------+----------------+-----------------+
| 29    |              | Pca21          | P c a 21       | `SG['Pca21']`   |
+-------+              +----------------+----------------+-----------------+
| 30    |              | Pnc2           | P n c 2        | `SG['Pnc2']`    |
+-------+              +----------------+----------------+-----------------+
| 31    |              | Pmn21          | P m n 21       | `SG['Pmn21']`   |
+-------+              +----------------+----------------+-----------------+
| 32    |              | Pba2           | P b a 2        | `SG['Pba2']`    |
+-------+              +----------------+----------------+-----------------+
| 33    |              | Pna21          | P n a 21       | `SG['Pna21']`   |
+-------+              +----------------+----------------+-----------------+
| 34    |              | Pnn2           | P n n 2        | `SG['Pnn2']`    |
+-------+              +----------------+----------------+-----------------+
| 35    |              | Cmm2           | C m m 2        | `SG['Cmm2']`    |
+-------+              +----------------+----------------+-----------------+
| 36    |              | Cmc21          | C m c 21       | `SG['Cmc21']`   |
+-------+              +----------------+----------------+-----------------+
| 37    |              | Ccc2           | C c c 2        | `SG['Ccc2']`    |
+-------+              +----------------+----------------+-----------------+
| 38    |              | Amm2           | A m m 2        | `SG['Amm2']`    |
+-------+              +----------------+----------------+-----------------+
| 39    |              | Aem2           | A e m 2        | `SG['Aem2']`    |
+-------+              +----------------+----------------+-----------------+
| 40    |              | Ama2           | A m a 2        | `SG['Ama2']`    |
+-------+              +----------------+----------------+-----------------+
| 41    |              | Aea2           | A e a 2        | `SG['Aea2']`    |
+-------+              +----------------+----------------+-----------------+
| 42    |              | Fmm2           | F m m 2        | `SG['Fmm2']`    |
+-------+              +----------------+----------------+-----------------+
| 43    |              | Fdd2           | F d d 2        | `SG['Fdd2']`    |
+-------+              +----------------+----------------+-----------------+
| 44    |              | Imm2           | I m m 2        | `SG['Imm2']`    |
+-------+              +----------------+----------------+-----------------+
| 45    |              | Iba2           | I b a 2        | `SG['Iba2']`    |
+-------+              +----------------+----------------+-----------------+
| 46    |              | Ima2           | I m a 2        | `SG['Ima2']`    |
+-------+              +----------------+----------------+-----------------+
| 47-74 |              | To             | be             | defined         |
+-------+--------------+----------------+----------------+-----------------+
| 75    | tetragonal   | To             | be             | defined         |
+-------+              +----------------+----------------+-----------------+
| -142  |              | To             | be             | defined         |
+-------+--------------+----------------+----------------+-----------------+
| 143-  | trigonal     | To             | be             | defined         |
+-------+              +----------------+----------------+-----------------+
| -167  |              | To             | be             | defined         |
+-------+--------------+----------------+----------------+-----------------+
| 168-  | hexagonal    | To             | be             | defined         |
+-------+              +----------------+----------------+-----------------+
| -194  |              | To             | be             | defined         |
+-------+--------------+----------------+----------------+-----------------+
| 195-  | cubic        | To             | be             | defined         |
+-------+              +----------------+----------------+-----------------+
| -230  |              | To             | be             | defined         |
+-------+--------------+----------------+----------------+-----------------+
"""

if __name__ == '__main__':
    g = SG['Aea2']; print(g); [print(op) for op in g]
    g = SG['Fmm2']; print(g); [print(op) for op in g]
    g = SG['Fdd2']; print(g); [print(op) for op in g]
    g = SG['Imm2']; print(g); [print(op) for op in g]
    g = SG['Iba2']; print(g); [print(op) for op in g]
    g = SG['Ima2']; print(g); [print(op) for op in g]
