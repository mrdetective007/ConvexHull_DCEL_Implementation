/** @file Face.h
 * @brief This file contains the definition of the Face struct.
 * @details The Face struct represents a face in a planar subdivision.
 *
 */
#ifndef FACE_H
#define FACE_H

#include <stddef.h>
#include "HalfEdge.h"
using namespace std;

struct HalfEdge;

// This struct represents a face in a planar subdivision
struct Face
{
public:
    /** @class Face
     * @brief This struct represents a face in a planar subdivision.
     * It stores a pointer to a half-edge that lies on the boundary of the face and a boolean flag indicating whether the face is still valid or not.
     * @param edge Pointer to a half-edge that lies on the boundary of the face
     * @param is_valid A boolean flag indicating whether the face is still valid or not
     * @note The half-edges are stored in a doubly connected edge list (DCEL) data structure.
     * @note The face is valid if it is not intersected by any of the edges of the subdivision.
     * @note The face is invalid if it is intersected by any of the edges of the subdivision.
     */

    // Pointer to a half-edge that lies on the boundary of the face
    HalfEdge *edge;

    // A boolean flag indicating whether the face is still valid or not
    bool is_valid;
    Face();
};

#endif // FACE_H
