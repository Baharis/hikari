import numpy as np
import random


def fibonacci_sphere(samples=1, seed=1337):
    random.seed(seed)
    rnd = random.random() * samples
    points = []
    offset = 2. / samples
    increment = np.pi * (3. - np.sqrt(5.))
    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - pow(y, 2))
        phi = ((i + rnd) % samples) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append((x, y, z))
    return points
