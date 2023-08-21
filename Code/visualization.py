"""
Visualization of the input polygon, output polygon and merged polygon
"""
import matplotlib.pyplot as plt
import numpy as np
 
fig, ax = plt.subplots(1,3,sharey=True)
 
f = open('input.txt', 'r')
polygon = []
n = int((f.readline()))
for i in range(n):
    polygon.append([float(x) for x in f.readline().split()])
f.close()
x = []
y = []
for i in range(len(polygon)):
    x.append(polygon[i][0])
    y.append(polygon[i][1])
x.append(polygon[0][0])
y.append(polygon[0][1])
ax[0].set_title("Number of faces in the input polygon : "+ " 1")
ax[0].plot(x, y, 'r-')
 
f = open('output.txt', 'r')
polygons = []
n = int(f.readline())
for i in range(n):
    polygon = []
    m = int(f.readline())
    for j in range(m):
        polygon.append([float(x) for x in f.readline().split()])
    polygons.append(polygon)
f.close()
for polygon in polygons:
    x = []
    y = []
    for i in range(len(polygon)):
        x.append(polygon[i][0])
        y.append(polygon[i][1])
    x.append(polygon[0][0])
    y.append(polygon[0][1])
    ax[1].set_title("Number of faces in the output polygon : " + str(n))
    ax[1].plot(x, y)
 
 
f = open('merge.txt', 'r')
polygons = []
nn = int(f.readline())
for i in range(nn):
    polygon = []
    mm = int(f.readline())
    for j in range(mm):
        polygon.append([float(x) for x in f.readline().split()])
    polygons.append(polygon)
f.close()
for polygon in polygons:
    x = []
    y = []
    for i in range(len(polygon)):
        x.append(polygon[i][0])
        y.append(polygon[i][1])
    x.append(polygon[0][0])
    y.append(polygon[0][1])
    ax[2].set_title("Number of faces in the merged polygon :" + str(nn))
    ax[2].plot(x, y)
 
plt.show()
