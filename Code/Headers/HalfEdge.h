/** @file Halfedge.h
 * @brief This file contains the definition of the HalfEdge struct.
 * @details The HalfEdge struct represents a half-edge in a planar subdivision.
 */
#ifndef HALFEDGE_H
#define HALFEDGE_H

#include <stddef.h>
#include "Vertex.h" // assuming Vertex class is defined in vertex.h
#include "Face.h"   // assuming Face class is defined in face.h
using namespace std;

struct Face;
struct Vertex;

struct HalfEdge
{
public:
    /**
     * @class HalfEdge
     * @brief This struct represents a half-edge in a planar subdivision.
     * It stores pointers to the origin vertex, the next half-edge in the face, the previous half-edge in the face, the twin half-edge and the face the half-edge belongs to.
     * @param origin Pointer to the origin vertex of the half-edge
     * @param next Pointer to the next half-edge in the face
     * @param prev Pointer to the previous half-edge in the face
     * @param twin Pointer to the twin half-edge (opposite direction)
     * @param face Pointer to the face the half-edge belongs to
     * @note The half-edges are stored in a doubly connected edge list (DCEL) data structure.
     */

    Vertex *origin; // pointer to the origin vertex of the half-edge
    HalfEdge *next; // pointer to the next half-edge in the face
    HalfEdge *prev; // pointer to the previous half-edge in the face
    HalfEdge *twin; // pointer to the twin half-edge (opposite direction)
    Face *face;     // pointer to the face the half-edge belongs to

    /**
     * @brief Constructor for HalfEdge struct
     * @param origin Pointer to the origin vertex of the half-edge
     * @param face Pointer to the face the half-edge belongs to
     */
    HalfEdge(Vertex *origin, Face *face);
};

#endif /* HALFEDGE_H */
