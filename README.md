# Convex Hull and DCEL Implementation

Welcome to the Convex Hull and DCEL Implementation repository! This repository contains C++ code files for implementing the Convex Hull algorithm and the Doubly Connected Edge List (DCEL) data structure.

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Code Files](#code-files)
- [Usage](#usage)
- [Explanation of Code Files](#explanation-of-code-files)
  - [ConvexHull.cpp](#convexhullcpp)
  - [DCEL.cpp](#dcelcpp)
- [Conclusion](#conclusion)

## Introduction

This repository showcases the implementation of two fundamental concepts in computational geometry: the Convex Hull algorithm and the Doubly Connected Edge List (DCEL) data structure. The Convex Hull algorithm is used to find the convex hull of a set of points in a plane, while the DCEL data structure represents planar subdivisions.

## Prerequisites

To compile and run the C++ code files, you'll need:

- C++ compiler (e.g., g++)
- Basic knowledge of computational geometry and data structures

## Code Files

This repository contains the following C++ code files:

- [ConvexHull.cpp](Code/ConvexHull.cpp): Implementation of the Convex Hull algorithm using the Graham's Scan technique.
- [DCEL.cpp](Code/DCEL.cpp): Implementation of the Doubly Connected Edge List (DCEL) data structure.

## Usage

1. Clone the repository:

```bash
git clone https://github.com/mrdetective007/ConvexHull_DCEL_Implementation.git
```

2. Navigate to the `Code` directory:

```bash
cd ConvexHull_DCEL_Implementation/Code
```

3. Compile and run the C++ code:

For Convex Hull:
```bash
g++ ConvexHull.cpp -o convex_hull
./convex_hull
```

For DCEL:
```bash
g++ DCEL.cpp -o dcel
./dcel
```

4. Follow the instructions in each code file to input the necessary data. The code will output the results of the Convex Hull algorithm or demonstrate the functionality of the DCEL.

## Explanation of Code Files

### ConvexHull.cpp

The `ConvexHull.cpp` file implements the Convex Hull algorithm using the Graham's Scan technique. The algorithm proceeds as follows:

1. Input a set of points in the plane.
2. Find the point with the lowest y-coordinate (and the leftmost if tied) as the pivot.
3. Sort the remaining points based on their polar angles with respect to the pivot.
4. Initialize an empty stack to hold the convex hull points.
5. Iterate through the sorted points and perform a series of left or right turns to build the convex hull.

The code outputs the vertices of the convex hull in counterclockwise order.

### DCEL.cpp

The `DCEL.cpp` file implements the Doubly Connected Edge List (DCEL) data structure. The DCEL represents planar subdivisions, such as the result of the Convex Hull algorithm. It consists of several components:

1. Vertex: Represents a vertex in the planar subdivision.
2. Half-Edge: Represents an edge as two half-edges, each pointing in opposite directions.
3. Face: Represents a face in the planar subdivision.

The code demonstrates the functionalities of adding vertices, edges, and faces, as well as traversing the structure.

## Conclusion

The Convex Hull algorithm and the Doubly Connected Edge List (DCEL) data structure are foundational concepts in computational geometry. The provided code files showcase their implementations and applications in understanding geometric algorithms and data structures.
