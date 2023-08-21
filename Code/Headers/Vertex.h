
/** @file Vertex.h
 * @brief This file contains the definition of the Vertex struct.
 * @details The Vertex struct represents a vertex in a planar subdivision.
 */
#ifndef VERTEX_H
#define VERTEX_H

#include <stddef.h>
#include "HalfEdge.h" // assuming HalfEdge class is defined in Halfedge.h
using namespace std;

struct HalfEdge;

struct Vertex
{
public:
    /**
     * @class Vertex
     * @brief This struct represents a vertex in a planar subdivision.
     * It stores the (x,y) coordinates of the vertex and a pointer to an outgoing half-edge.
     * @note The outgoing half-edge is the half-edge that is incident on the vertex and lies on the inside of the polygon.
     * @param x The x-coordinate of the vertex
     * @param y The y-coordinate of the vertex
     * @param id The id of the vertex
     * @param edge Pointer to the outgoing half-edge
     *
     */

    double x;
    double y;
    HalfEdge *edge; // outgoing edge corresponding to the inside polygon
    int id;

    /**
     * @brief Constructor for Vertex struct
     * @param x The x-coordinate of the vertex
     * @param y The y-coordinate of the vertex
     * @param id The id of the vertex
     */
    Vertex(double x, double y, int id);
};

#endif /* VERTEX_H */
