/** @file DCEL.h
 * @brief This file contains the definition of the DCEL class.
 * @details The DCEL class represents a Doubly Connected Edge List (DCEL) data structure.
 */
#ifndef DCEL_H
#define DCEL_H

#include <vector>
#include "HalfEdge.h"
#include "Vertex.h"
#include "Face.h"

struct Face;
struct HalfEdge;
struct Vertex;

class DCEL
{
public:
    /** @class DCEL
     * @brief This class represents a Doubly Connected Edge List (DCEL) data structure.
     * It stores the vertices, faces and half-edges of a planar subdivision.
     * It also provides functions to create a polygon and traverse the faces of the subdivision.
     * @note The DCEL data structure is used to represent a planar subdivision.
     */
    std::vector<Vertex *> vertices; // vector to store vertices
    std::vector<Face *> faces;      // vector to store faces

    std::vector<HalfEdge *> diagonals; // vector to store diagonals
    std::vector<HalfEdge *> edges;     // vectors to store edges

    // Constructor
    DCEL();

    // Function to create a polygon given a set of vertices
    void create_polygon(std::vector<std::pair<double, double>> &arr);
};

#endif // DCEL_H