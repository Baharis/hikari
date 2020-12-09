from kesshou.symmetry import Group, SymmOp
import numpy as np
from pathlib import Path
import json

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
az = SymmOp.from_code('x+1/2, y, -z')
bx = SymmOp.from_code('-x, y+1/2, z')
cx = SymmOp.from_code('-x, y, z+1/2')
cy = SymmOp.from_code('x, -y, z+1/2')
dx = SymmOp.from_code('-x, y+1/4, z+1/4')
dy = SymmOp.from_code('x+1/4, -y, z+1/4')
dz = SymmOp.from_code('x+1/4, y+1/4, -z')
mx = SymmOp.from_code('-x, y, z')
my = SymmOp.from_code('x, -y, z')
mz = SymmOp.from_code('x, y, -z')
nx = SymmOp.from_code('-x, y+1/2, z+1/2')
ny = SymmOp.from_code('x+1/2, -y, z+1/2')
nz = SymmOp.from_code('x+1/2, y+1/2, -z')
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
    'Fdd2': Group(A, B, r2z, dx.at(x/2)),
    'Imm2': Group(I, r2z, mx),
    'Iba2': Group(I, r2z, cx),
    'Ima2': Group(I, r2z, ay),
    'Pmmm': Group(mx, my, mz),
    'Pnnn': Group(nx.at(x), ny.at(y), nz.at(z)),
    'Pccm': Group(cx, cy, mz),
    'Pban': Group(bx.at(x), ay.at(y), nz),
    'Pmma': Group(mx.at(x), my, az),
    'Pnna': Group(nx, ny.at(y), az),
    'Pmna': Group(mx, ny, az.at(z)),
    'Pcca': Group(cx.at(x), cy, az),
    'Pbam': Group(bx.at(x), ay.at(y), mz),
    'Pccn': Group(cx.at(x), cy.at(y), nz),
    'Pbcm': Group(bx, cy.at(y), mz.at(z)),
    'Pnnm': Group(nx.at(x), ny.at(y), mz),
    'Pmmn': Group(mx, my, nz),
    'Pbcn': Group(bx.at(x), cy, nz.at(z)),
    'Pbca': Group(bx.at(x), cy.at(y), az.at(z)),
    'Pnma': Group(nx.at(x), my.at(y), az.at(z)),
    'Cmcm': Group(C, mx, cy, mz.at(z)),
    'Cmca': Group(C, mx, cy.at(y), az.at(z)),
    'Cmmm': Group(C, mx, my, mz),
    'Cccm': Group(C, cx, cy, mz),
    'Cmme': Group(C, mx, my.at(y), az),
    'Ccce': Group(C, cx.at(x), cy.at(y), az.at(z)),
    'Fmmm': Group(A, B, mx, my, mz),
    'Fddd': Group(A, B, dx.at(x/2), dy.at(y/2), dz.at(z/2)),
    'Immm': Group(I, mx, my, mz),
    'Ibam': Group(I, bx.at(x), ay.at(x), mz),
    'Ibca': Group(I, bx.at(x), cy.at(y), az.at(z)),
    'Imma': Group(I, mx, my.at(y), az.at(z))
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

+-------+--------------+----------------+-------------------+-----------------+
| No.   | CRYSTAL      | Hermann-Maugin | Hermann-Maugin    | Can be accessed |
|       | SYSTEM       | short notation | full notation     | using           |
+-------+--------------+----------------+-------------------+-----------------+
| 1     | triclinic    | P1             | P 1               | `SG['P1']`      |
+-------+              +----------------+-------------------+-----------------+
| 2     |              | P-1            | P -1              | `SG['P-1']`     |
+-------+--------------+----------------+-------------------+-----------------+
| 3     | monoclinic   | P2             | P 1 2 1           | `SG['P2']`      |
+-------+              +----------------+-------------------+-----------------+
| 4     |              | P21            | P 1 21 1          | `SG['P21']`     |
+-------+              +----------------+-------------------+-----------------+
| 5     |              | C2             | C 1 2 1           | `SG['C2']`      |
+-------+              +----------------+-------------------+-----------------+
| 6     |              | Pm             | P 1 m 1           | `SG['Pm']`      |
+-------+              +----------------+-------------------+-----------------+
| 7     |              | Pc             | P 1 c 1           | `SG['Pc']`      |
+-------+              +----------------+-------------------+-----------------+
| 8     |              | Cm             | C 1 m 1           | `SG['Cm']`      |
+-------+              +----------------+-------------------+-----------------+
| 9     |              | Cc             | C 1 c 1           | `SG['Cc']`      |
+-------+              +----------------+-------------------+-----------------+
| 10    |              | P2/m           | P 1 2/m 1         | `SG['P2/m']`    |
+-------+              +----------------+-------------------+-----------------+
| 11    |              | P21/m          | P 1 21/m 1        | `SG['P21/m']`   |
+-------+              +----------------+-------------------+-----------------+
| 12    |              | C2/m           | C 1 2/m 1         | `SG['C2/m']`    |
+-------+              +----------------+-------------------+-----------------+
| 13    |              | P2/c           | P 1 2/c 1         | `SG['P2/c']`    |
+-------+              +----------------+-------------------+-----------------+
| 14    |              | P21/c          | P 1 21/c 1        | `SG['P21/c']`   |
+-------+              +----------------+-------------------+-----------------+
| 15    |              | C2/c           | C 1 2/c 1         | `SG['C2/c']`    | 
+-------+--------------+----------------+-------------------+-----------------+
| 16    | orthorhombic | P222           | P 2 2 2           | `SG['P222']`    |
+-------+              +----------------+-------------------+-----------------+
| 17    |              | P2221          | P 2 2 21          | `SG['P2221']`   |
+-------+              +----------------+-------------------+-----------------+
| 18    |              | P21212         | P 21 21 2         | `SG['P21212']`  |
+-------+              +----------------+-------------------+-----------------+
| 19    |              | P212121        | P 21 21 21        | `SG['P212121']` |
+-------+              +----------------+-------------------+-----------------+
| 20    |              | C2221          | C 2 2 21          | `SG['C2221']`   |
+-------+              +----------------+-------------------+-----------------+
| 21    |              | C222           | C 2 2 2           | `SG['C222']`    |
+-------+              +----------------+-------------------+-----------------+
| 22    |              | F222           | F 2 2 2           | `SG['F222']`    |
+-------+              +----------------+-------------------+-----------------+
| 23    |              | I222           | I 2 2 2           | `SG['I222']`    |
+-------+              +----------------+-------------------+-----------------+
| 24    |              | I212121        | I 21 21 21        | `SG['I212121']` |
+-------+              +----------------+-------------------+-----------------+
| 25    |              | Pmm2           | P m m 2           | `SG['Pmm2']`    |
+-------+              +----------------+-------------------+-----------------+
| 26    |              | Pmc21          | P m c 21          | `SG['Pmc21']`   |
+-------+              +----------------+-------------------+-----------------+
| 27    |              | Pcc2           | P c c 2           | `SG['Pcc2']`    |
+-------+              +----------------+-------------------+-----------------+
| 28    |              | Pma2           | P m a 2           | `SG['Pma2']`    |
+-------+              +----------------+-------------------+-----------------+
| 29    |              | Pca21          | P c a 21          | `SG['Pca21']`   |
+-------+              +----------------+-------------------+-----------------+
| 30    |              | Pnc2           | P n c 2           | `SG['Pnc2']`    |
+-------+              +----------------+-------------------+-----------------+
| 31    |              | Pmn21          | P m n 21          | `SG['Pmn21']`   |
+-------+              +----------------+-------------------+-----------------+
| 32    |              | Pba2           | P b a 2           | `SG['Pba2']`    |
+-------+              +----------------+-------------------+-----------------+
| 33    |              | Pna21          | P n a 21          | `SG['Pna21']`   |
+-------+              +----------------+-------------------+-----------------+
| 34    |              | Pnn2           | P n n 2           | `SG['Pnn2']`    |
+-------+              +----------------+-------------------+-----------------+
| 35    |              | Cmm2           | C m m 2           | `SG['Cmm2']`    |
+-------+              +----------------+-------------------+-----------------+
| 36    |              | Cmc21          | C m c 21          | `SG['Cmc21']`   |
+-------+              +----------------+-------------------+-----------------+
| 37    |              | Ccc2           | C c c 2           | `SG['Ccc2']`    |
+-------+              +----------------+-------------------+-----------------+
| 38    |              | Amm2           | A m m 2           | `SG['Amm2']`    |
+-------+              +----------------+-------------------+-----------------+
| 39    |              | Aem2           | A e m 2           | `SG['Aem2']`    |
+-------+              +----------------+-------------------+-----------------+
| 40    |              | Ama2           | A m a 2           | `SG['Ama2']`    |
+-------+              +----------------+-------------------+-----------------+
| 41    |              | Aea2           | A e a 2           | `SG['Aea2']`    |
+-------+              +----------------+-------------------+-----------------+
| 42    |              | Fmm2           | F m m 2           | `SG['Fmm2']`    |
+-------+              +----------------+-------------------+-----------------+
| 43    |              | Fdd2           | F d d 2           | `SG['Fdd2']`    |
+-------+              +----------------+-------------------+-----------------+
| 44    |              | Imm2           | I m m 2           | `SG['Imm2']`    |
+-------+              +----------------+-------------------+-----------------+
| 45    |              | Iba2           | I b a 2           | `SG['Iba2']`    |
+-------+              +----------------+-------------------+-----------------+
| 46    |              | Ima2           | I m a 2           | `SG['Ima2']`    |
+-------+              +----------------+-------------------+-----------------+
| 47    |              | Pmmm           | P 2/m 2/m 2/m     | `SG['Pmmm']`    |
+-------+              +----------------+-------------------+-----------------+
| 48    |              | Pnnn           | P 2/n 2/n 2/n     | `SG['Pnnn']`    |
+-------+              +----------------+-------------------+-----------------+
| 49    |              | Pccm           | P 2/c 2/c 2/m     | `SG['Pccm']`    |
+-------+              +----------------+-------------------+-----------------+
| 50    |              | Pban           | P 2/b 2/a 2/n     | `SG['Pban']`    |
+-------+              +----------------+-------------------+-----------------+
| 51    |              | Pmma           | P 21/m 2/m 2/a    | `SG['Pmma']`    |
+-------+              +----------------+-------------------+-----------------+
| 52    |              | Pnna           | P 2/n 21/n 2/a    | `SG['Pnna']`    |
+-------+              +----------------+-------------------+-----------------+
| 53    |              | Pmna           | P 2/m 2/n 21/a    | `SG['Pmna']`    |
+-------+              +----------------+-------------------+-----------------+
| 54    |              | Pcca           | P 21/c 2/c 2/a    | `SG['Pcca']`    |
+-------+              +----------------+-------------------+-----------------+
| 55    |              | Pbam           | P 21/b 21/a 2/m   | `SG['Pbam']`    |
+-------+              +----------------+-------------------+-----------------+
| 56    |              | Pccn           | P 21/c 21/c 2/n   | `SG['Pccn']`    |
+-------+              +----------------+-------------------+-----------------+
| 57    |              | Pbcm           | P 2/b 21/c 21/m   | `SG['Pbcm']`    |
+-------+              +----------------+-------------------+-----------------+
| 58    |              | Pnnm           | P 21/n 21/n 2/m   | `SG['Pnnm']`    |
+-------+              +----------------+-------------------+-----------------+
| 59    |              | Pmmn           | P 21/m 21/m 2/n   | `SG['Pmmn']`    |
+-------+              +----------------+-------------------+-----------------+
| 60    |              | Pbcn           | P 21/b 2/c 21/n   | `SG['Pbcn']`    |
+-------+              +----------------+-------------------+-----------------+
| 61    |              | Pbca           | P 21/b 21/c 21/n  | `SG['Pbca']`    |
+-------+              +----------------+-------------------+-----------------+
| 62    |              | Pnma           | P 21/n 21/m 21/a  | `SG['Pnma']`    |
+-------+              +----------------+-------------------+-----------------+
| 63    |              | Cmcm           | C 2/m 2/c 21/m    | `SG['Cmcm']`    |
+-------+              +----------------+-------------------+-----------------+
| 64    |              | Cmca           | C 2/m 2/c 21/a    | `SG['Cmca']`    |
+-------+              +----------------+-------------------+-----------------+
| 65    |              | Cmmm           | C 2/m 2/m 2/m     | `SG['Cmmm']`    |
+-------+              +----------------+-------------------+-----------------+
| 66    |              | Cccm           | C 2/c 2/c 2/m     | `SG['Cccm']`    |
+-------+              +----------------+-------------------+-----------------+
| 67    |              | Cmme           | C 2/m 2/m 2/e     | `SG['Cmme']`    |
+-------+              +----------------+-------------------+-----------------+
| 68    |              | Ccce           | C 2/c 2/c 2/e     | `SG['Ccce']`    |
+-------+              +----------------+-------------------+-----------------+
| 69    |              | Fmmm           | F 2/m 2/m 2/m     | `SG['Fmmm']`    |
+-------+              +----------------+-------------------+-----------------+
| 70    |              | Fddd           | F 2/d 2/d 2/d     | `SG['Fddd']`    |
+-------+              +----------------+-------------------+-----------------+
| 71    |              | Immm           | I 2/m 2/m 2/m     | `SG['Immm']`    |
+-------+              +----------------+-------------------+-----------------+
| 72    |              | Ibam           | I 2/b 2/a 2/m     | `SG['Ibam']`    |
+-------+              +----------------+-------------------+-----------------+
| 73    |              | Ibca           | I 2/b 2/c 2/a     | `SG['Ibca']`    |
+-------+              +----------------+-------------------+-----------------+
| 74    |              | Imma           | I 2/m 2/m 2/a     | `SG['Imma']`    |
+-------+--------------+----------------+-------------------+-----------------+
| 75    | tetragonal   | To             | be                | defined         |
+-------+              +----------------+-------------------+-----------------+
| -142  |              | To             | be                | defined         |
+-------+--------------+----------------+-------------------+-----------------+
| 143-  | trigonal     | To             | be                | defined         |
+-------+              +----------------+-------------------+-----------------+
| -167  |              | To             | be                | defined         |
+-------+--------------+----------------+-------------------+-----------------+
| 168-  | hexagonal    | To             | be                | defined         |
+-------+              +----------------+-------------------+-----------------+
| -194  |              | To             | be                | defined         |
+-------+--------------+----------------+-------------------+-----------------+
| 195-  | cubic        | To             | be                | defined         |
+-------+              +----------------+-------------------+-----------------+
| -230  |              | To             | be                | defined         |
+-------+--------------+----------------+-------------------+-----------------+
"""


def _unpack_space_group_dictionary_from_json():
    current_file_path = Path(__file__).parent.absolute()
    json_file_path = current_file_path.joinpath('space_groups.json')
    with open(json_file_path) as file:
        json_dict = json.load(file)
    space_group_dict = {}
    for _, json_group in json_dict.items():
        new_group_name = json_group["H-M_short"]
        print(json_group["generators"])
        new_group_gens = [SymmOp.from_code(c) for c in json_group["generators"]]
        new_group_object = Group(*new_group_gens)
        space_group_dict[new_group_name] = new_group_object
    return space_group_dict


SG2 = _unpack_space_group_dictionary_from_json()

if __name__ == '__main__':
    print(SG2)
    g = SG2['I41']; print(g); [print(op) for op in g]
    g = SG2['P-4']; print(g); [print(op) for op in g]
    g = SG2['I-1']; print(g); [print(op) for op in g]
    g = SG2['P4/m']; print(g); [print(op) for op in g]

