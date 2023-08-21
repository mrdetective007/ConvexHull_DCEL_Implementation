"""
Generate a random polygon with a given number of vertices.
"""
import numpy as np
import matplotlib.pyplot as plt
import random
from polygenerator import (
    random_polygon,
)


def plot_polygon(polygon):
    plt.figure()
    plt.gca().set_aspect("equal")
    polygon.append(polygon[0])
    xs, ys = zip(*polygon)
    plt.plot(xs, ys, "b-", linewidth=1)
    plt.show()


needseed = random.randint(0, 1000)
random.seed(needseed)
# num_points = random.randint(10, 100)
num_points = 4267
polygon = random_polygon(num_points)
polygon.reverse()
# plot_polygon(polygon)
print(num_points)
for i in range(len(polygon)):
    print(round(polygon[i][0], 5), round(polygon[i][1], 5))
