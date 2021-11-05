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
|       |               | -42m           | D2d           | `PG['-42m']`    |
| 14    |               +----------------+---------------+-----------------+
|       |               | -4m2           | D2d           | `PG['-4m2']`    |
+-------+               +----------------+---------------+-----------------+
| 15    |               | 4/mmm          | D4h           | `PG['4/mmm']`   |
+-------+---------------+----------------+---------------+-----------------+
| 16    | trigonal      | 3              | C3            | `PG['3']`       |
+-------+               +----------------+---------------+-----------------+
| 17    |               | -3             | C3i           | `PG['-3']`      |
+-------+               +----------------+---------------+-----------------+
|       |               | 321            | D3            | `PG['321']`     |
| 18    |               +----------------+---------------+-----------------+
|       |               | 312            | D3            | `PG['312']`     |
+-------+               +----------------+---------------+-----------------+
|       |               | 3m1            | C3v           | `PG['3m1']`     |
| 19    |               +----------------+---------------+-----------------+
|       |               | 31m            | C3v           | `PG['31m']`     |
+-------+               +----------------+---------------+-----------------+
|       |               | -3m1           | D3d           | `PG['-3m1']`    |
| 20    |               +----------------+---------------+-----------------+
|       |               | -31m           | D3d           | `PG['-31m']`    |
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
|       |               | -6m2           | D3h           | `PG['-6m2']`    |
| 26    |               +----------------+---------------+-----------------+
|       |               | -62m           | D3h           | `PG['-62m']`    |
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
"""
