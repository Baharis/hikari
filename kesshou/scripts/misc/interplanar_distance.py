import numpy as np
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import functools


# zdefiniuj zmienne
molecule1_file_path = '/home/dtchon/x/AP/interplanar_distance/21_mol5.tsv'
molecule2_file_path = '/home/dtchon/x/AP/interplanar_distance/21_mol6.tsv'

# w plikach powinny być jedynie pozycje atomów z jednej i drugiej cząsteczki
# wczytaj i zinterpretuj pliki
molecule1 = np.loadtxt(molecule1_file_path)
molecule2 = np.loadtxt(molecule2_file_path)

# znajdź centrum każdej cząsteczki i odległość między nimi
centre1 = np.mean(molecule1, axis=0)
centre2 = np.mean(molecule2, axis=0)


# przesuń cząsteczkę 1 w pozycję drugiej
molecule1 = molecule1 - centre1
molecule2 = molecule2 - centre2


# zdefiniuj płaszczyzną do dofittowania
def plane(x, y, pars):
    a = pars[0]
    b = pars[1]
    c = pars[2]
    z = a*x + b*y + c          # 0 = ax + by -1z + c
    return z


# zdefiniuj funkcję błędu
def error(pars, points):
    result = 0
    for x, y, z in points:
        plane_z = plane(x, y, pars)
        diff = abs(plane_z - z)
        result += diff**2
    return result


# zdefiniuj trójwymiarowy produkt wektorów
def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]


# przygotuj wspólną cząsteczkę, funkcję błędu, punkt początkowy
molecules = np.concatenate((molecule1, molecule2), axis=0)
fun = functools.partial(error, points=molecules)
pars0 = np.array([0, 0, 0])

# dofittuj płaszczyznę i zapisz wynik
res = opt.minimize(fun, pars0)
a, b, c = res.x[0:3]

# narysuj punkty w przestrzeni
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(*zip(*molecules))

# wektor normalny do płaszczyzny
point = np.array([0.0, 0.0, c])
normal = np.array(cross([1, 0, a], [0, 1, b]))

# odległość między centroidami to
vector = centre2 - centre1
vector_magnitude = np.linalg.norm(vector)
print('Odległość między centroidami:   ', vector_magnitude, sep='\t')

# odległość między płaszczyznami to
interplanar = np.dot(normal / np.linalg.norm(normal), vector)
interplanar_magnitude = np.linalg.norm(interplanar)
print('Odległość między płaszczyznami: ', interplanar_magnitude, sep='\t')

# odległość horyzontalna centroid to
horizontal_magnitude = np.sqrt(vector_magnitude**2 - interplanar_magnitude**2)
print('Odległość horyzontalna centroid:', horizontal_magnitude, sep='\t')

# TODO add documentation and translate to english
