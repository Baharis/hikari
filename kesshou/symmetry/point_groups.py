from kesshou.symmetry import Group, SymmOp

PG = {
    # TRICLINIC
    '1': Group(SymmOp.from_code('x,y,z')),
    '-1': Group(SymmOp.from_code('-x,-y,-z')),
    # MONOCLINIC
    '2': Group(SymmOp.from_code('-x,y,-z')),
    'm': Group(SymmOp.from_code('x,-y,z')),
    '2/m': Group(SymmOp.from_code('-x,-y,-z'), SymmOp.from_code('x,-y,z')),
    # ORTHORHOMBIC
    '222': Group(SymmOp.from_code('x,-y,-z'), SymmOp.from_code('-x,y,-z')),
    'mm2': Group(SymmOp.from_code('-x,y,z'), SymmOp.from_code('-x,-y,z')),
    'mmm': Group(SymmOp.from_code('-x,y,z'), SymmOp.from_code('x,-y,z'),
                 SymmOp.from_code('x,y,-z')),
    # TETRAGONAL
    '4': Group(SymmOp.from_code('-y,x,z')),
    '-4': Group(SymmOp.from_code('y,-x,-z')),
    '4/m': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('x,y,-z')),
    '422': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('-x,y,-z')),
    '4mm': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('x,-y,z')),
    '-42m': Group(SymmOp.from_code('y,-x,-z'), SymmOp.from_code('x,-y,-z')),
    '-4m2': Group(SymmOp.from_code('y,-x,-z'), SymmOp.from_code('-x,y,z')),
    '4/mmm': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('x,y,-z'),
                   SymmOp.from_code('-x,y,z')),
    # TRIGONAL
    '3': Group(SymmOp.from_code('-y,x-y,z')),
    '-3': Group(SymmOp.from_code('y,-x+y,-z')),
    '321': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('y,x,-z')),
    '312': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('-y,-x,-z')),
    '3m1': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('-y,-x,z')),
    '31m': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('y,x,z')),
    '-3m1': Group(SymmOp.from_code('y,-x+y,-z'), SymmOp.from_code('-y,-x,z')),
    '-31m': Group(SymmOp.from_code('y,-x+y,-z'), SymmOp.from_code('y,x,z')),
    # HEXAGONAL
    '6': Group(SymmOp.from_code('x-y,x,z')),
    '-6': Group(SymmOp.from_code('-x+y,-x,-z')),
    '6/m': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('x,y,-z')),
    '622': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('x-y,-y,-z')),
    '6mm': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('-x+y,y,z')),
    '-6m2': Group(SymmOp.from_code('-x+y,-x,-z'), SymmOp.from_code('-x+y,y,z')),
    '-62m': Group(SymmOp.from_code('-x+y,-x,-z'), SymmOp.from_code('y,x,-z')),
    '6/mmm': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('x,y,-z'),
                   SymmOp.from_code('-x+y,y,z')),
    # CUBIC
    '23': Group(SymmOp.from_code('-x,-y,z'), SymmOp.from_code('z,x,y')),
    'm-3': Group(SymmOp.from_code('x,y,-z'), SymmOp.from_code('-z,-x,-y')),
    '432': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('z,x,y')),
    '-43m': Group(SymmOp.from_code('y,-x,-z'), SymmOp.from_code('z,x,y')),
    'm-3m': Group(SymmOp.from_code('x,y,-z'), SymmOp.from_code('-z,-x,-y'),
                  SymmOp.from_code('-y,-x,z'))}
"""
Dictionary containing all known point groups written as :class:`Group`
along with alternative axis settings. The point groups in this dictionary
can be accessed using their short Hermann-Maugin notation, as presented below.

+-------+---------------+----------------+---------------+-----------------+
| No.   | CRYSTAL       | Hermann-Maugin | Schoenflies   | Can be accessed |
|       | SYSTEM        | notation       | notation      | using           |
+-------+---------------+----------------+---------------+-----------------+
| 1     | triclinic     | 1              | C1            | `PG['1']`       |
+-------+               +----------------+---------------+-----------------+
| 2     |               | -1             | Ci            | `PG['-1']`      |
+-------+---------------+----------------+---------------+-----------------+
| 3     | monoclinic    | 2              | C2            | `PG['2']`       |
+-------+               +----------------+---------------+-----------------+
| 4     |               | m              | Cs            | `PG['m']`       |
+-------+               +----------------+---------------+-----------------+
| 5     |               | 2/m            | C2h           | `PG['2/m']`     |
+-------+---------------+----------------+---------------+-----------------+
| 6     | orthorhombic  | 222            | D2            | `PG['222']`     |
+-------+               +----------------+---------------+-----------------+
| 7     |               | mm2            | C2v           | `PG['mm2']`     |
+-------+               +----------------+---------------+-----------------+
| 8     |               | mmm            | D2h           | `PG['mmm']`     |
+-------+---------------+----------------+---------------+-----------------+
| 9     | tetragonal    | 4              | C4            | `PG['4']`       |
+-------+               +----------------+---------------+-----------------+
| 10    |               | -4             | S4            | `PG['-4']`      |
+-------+               +----------------+---------------+-----------------+
| 11    |               | 4/m            | C4h           | `PG['4/m']`     |
+-------+               +----------------+---------------+-----------------+
| 12    |               | 422            | D4            | `PG['422']`     |
+-------+               +----------------+---------------+-----------------+
| 13    |               | 4mm            | C4v           | `PG['4mm']`     |
+-------+               +----------------+---------------+-----------------+
| 14    |               | -42m           | D2d           | `PG['-42m']`    |
+-------+               +----------------+---------------+-----------------+
| 14*   |               | -4m2           | D2d           | `PG['-4m2']`    |
+-------+               +----------------+---------------+-----------------+
| 15    |               | 4/mmm          | D4h           | `PG['4/mmm']`   |
+-------+---------------+----------------+---------------+-----------------+
| 16    | trigonal      | 3              | C3            | `PG['3']`       |
+-------+               +----------------+---------------+-----------------+
| 17    |               | -3             | C3i           | `PG['-3']`      |
+-------+               +----------------+---------------+-----------------+
| 18    |               | 32             | D3            | `PG['321']`     |
+-------+               +----------------+---------------+-----------------+
| 18*   |               | 32             | D3            | `PG['312']`     |
+-------+               +----------------+---------------+-----------------+
| 19    |               | 3m             | C3v           | `PG['3m1']`     |
+-------+               +----------------+---------------+-----------------+
| 19*   |               | 3m             | C3v           | `PG['31m']`     |
+-------+               +----------------+---------------+-----------------+
| 20    |               | -3m            | D3d           | `PG['-3m1']`    |
+-------+               +----------------+---------------+-----------------+
| 20*   |               | -3m            | D3d           | `PG['-31m']`    |
+-------+---------------+----------------+---------------+-----------------+
| 21    | hexagonal     | 6              | C6            | `PG['6']`       |
+-------+               +----------------+---------------+-----------------+
| 22    |               | -6             | C3h           | `PG['-6']`      |
+-------+               +----------------+---------------+-----------------+
| 23    |               | 6/m            | C6h           | `PG['6/m']`     |
+-------+               +----------------+---------------+-----------------+
| 24    |               | 622            | D6            | `PG['622']`     |
+-------+               +----------------+---------------+-----------------+
| 25    |               | 6mm            | C6v           | `PG['6mm']`     |
+-------+               +----------------+---------------+-----------------+
| 26    |               | -6m2           | D3h           | `PG['-6m2']`    |
+-------+               +----------------+---------------+-----------------+
| 26*   |               | -62m           | D3h           | `PG['-62m']`    |
+-------+               +----------------+---------------+-----------------+
| 27    |               | 6/mmm          | D6h           | `PG['6/mmm']`   |
+-------+---------------+----------------+---------------+-----------------+
| 28    | cubic         | 23             | T             | `PG['23']`      |
+-------+               +----------------+---------------+-----------------+
| 29    |               | m-3            | Th            | `PG['m-3']`     |
+-------+               +----------------+---------------+-----------------+
| 30    |               | 432            | O             | `PG['432']`     |
+-------+               +----------------+---------------+-----------------+
| 31    |               | -43m           | Td            | `PG['-43m']`    |
+-------+               +----------------+---------------+-----------------+
| 32    |               | m-3m           | Oh            | `PG['m-3m']`    |
+-------+---------------+----------------+---------------+-----------------+

Asterisk (*) denotes alternative choice of axes.
"""