/** @file main.cpp
 * @brief This file contains the implementation of the Algorithm 1 of the paper "Algorithms for the decomposition of a polygon into convex polygons" by J. Fernandez, L. Canovas, B. Pelegrin
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <ctime>
#include <iomanip>

#include "Headers/Vertex.h"
#include "Headers/Face.h"
#include "Headers/DCEL.h"
#include "Headers/HalfEdge.h"

using namespace std;

struct Vertex;
struct Face;
struct HalfEdge;

Vertex::Vertex(double x, double y, int id)
{
    this->x = x;
    this->y = y;
    this->edge = NULL;
    this->id = id;
}

Face::Face()
{
    edge = nullptr;  // The face has no boundary edge initially
    is_valid = true; // The face is initially valid
}

HalfEdge::HalfEdge(Vertex *v, Face *f)
{
    this->origin = v;
    this->face = f;
    this->twin = nullptr;
    this->next = nullptr;
    this->prev = nullptr;
}

DCEL::DCEL()
{
    vertices.clear();
    faces.clear();
    edges.clear();
}

void DCEL::create_polygon(vector<pair<double, double>> &arr)
{
    /**
     * @brief This function creates a polygon given a set of vertices.
     * It creates a face and a half-edge for each vertex.
     * It then connects all the half-edges to form a doubly connected edge list.
     * @param arr A vector of pairs of doubles representing the (x,y) coordinates of the vertices of the polygon.
     * @return void
     */
    int n = arr.size();

    // Create face
    Face *f = new Face();
    faces.push_back(f);

    // Create vertices
    for (int i = 0; i < n; i++)
    {
        Vertex *v = new Vertex(arr[i].first, arr[i].second, i + 1);
        vertices.push_back(v);
    }

    // Create half edges
    for (int i = 0; i < n; i++)
    {
        HalfEdge *e = new HalfEdge(vertices[i], f);
        vertices[i]->edge = e;
        edges.push_back(e);
    }

    // Connect all half edges
    for (int i = 0; i < n; i++)
    {
        edges[i]->prev = edges[(i - 1 + n) % n];
        edges[i]->next = edges[(i + 1) % n];
    }

    // Initialize Face
    f->edge = edges[0];
}

void traverse_face(Face *f)
{
    /** @brief This function traverses the boundary of the given face 'f' in counter-clockwise order.
     * It prints the number of edges (half-edges) in the boundary of the face.
     * Then, it traverses the half-edges in the boundary and prints the (x,y) coordinates of their origins.
     * The traversal is performed using a do-while loop to ensure that at least one iteration is performed even if the face has no edges.
     @param f Pointer to the face to be traversed
     @return void
     */
    // Start at an arbitrary half-edge in the boundary of the face
    HalfEdge *e = f->edge;
    int count = 0;
    // Traverse the boundary and count the number of half-edges
    do
    {
        count++;
        e = e->next;
    } while (e != f->edge);

    // Print the number of edges in the boundary
    cout << count << endl;

    // Traverse the boundary again and print the coordinates of the origin of each half-edge
    do
    {
        cout << e->origin->x << " " << e->origin->y << endl;
        e = e->next;
    } while (e != f->edge);
}

void traverse_face_merge(Face *f, FILE *fp)
{
    /** @brief This function traverses the boundary of the given face 'f' in counter-clockwise order.
     * It prints the number of edges (half-edges) in the boundary of the face.
     * Then, it traverses the half-edges in the boundary and prints the (x,y) coordinates of their origins.
     * The traversal is performed using a do-while loop to ensure that at least one iteration is performed even if the face has no edges.
     * @param f Pointer to the face to be traversed
     * @param fp Pointer to the file to write the output
     * @return void
     * @note This function is used to merge the output of the traverse_face function for all faces into a single file.
     */
    // Start at an arbitrary half-edge in the boundary of the face
    HalfEdge *e = f->edge;
    int count = 0;
    // Traverse the boundary and count the number of half-edges
    do
    {
        count++;
        e = e->next;
    } while (e != f->edge);

    // Print the number of edges in the boundary
    fprintf(fp, "%d\n", count);

    // Traverse the boundary again and print the coordinates of the origin of each half-edge
    do
    {
        fprintf(fp, "%lf %lf\n", e->origin->x, e->origin->y);
        e = e->next;
    } while (e != f->edge);
}

// This function splits a face to add a new diagonal edge
void split_face(Vertex *start, Vertex *end, DCEL *polygon)
{
    /** @brief This function splits a face to add a new diagonal edge.
     * It creates two new faces and two new half-edges.
     * It then connects the two new half-edges as twins.
     * It then removes the old face and adds the two new faces to the faces vector of the DCEL object.
     * It then connects the two new half-edges to the two new faces.
     * @param start Pointer to the vertex at the start of the diagonal
     * @param end Pointer to the vertex at the end of the diagonal
     * @param polygon Pointer to the DCEL object representing the polygon
     * @return void
     */
    // Create two new faces
    Face *f1 = new Face();
    Face *f2 = new Face();

    // Create two new half edges
    HalfEdge *h1 = new HalfEdge(end, f1);
    HalfEdge *h2 = new HalfEdge(start, f2);

    // Connect the two new half edges as twins
    h1->twin = h2;
    h2->twin = h1;

    // Add the new diagonal to the diagonals vector of the DCEL object
    polygon->diagonals.push_back(h1);

    // Set the edges of the two new faces
    f2->edge = h2;
    f1->edge = h1;

    // Remove the old face from the faces vector of the DCEL object and add the two new faces
    polygon->faces.pop_back();
    polygon->faces.push_back(f1);
    polygon->faces.push_back(f2);

    // Connect the new diagonal h1 to the existing polygon
    h1->next = start->edge;
    HalfEdge *temp = start->edge;
    while (temp->next->origin != end)
    {
        temp = temp->next;
    }
    h1->prev = temp;
    temp->next = h1;
    start->edge->prev = h1;

    temp = start->edge;
    // Allocate the new face to each edge of the new polygon
    while (temp != h1)
    {
        temp->face = f1;
        temp = temp->next;
    }

    // Connect the new diagonal h2 to the existing polygon
    h2->next = end->edge;
    temp = end->edge;
    while (temp->next->origin != start)
    {
        temp = temp->next;
    }
    h2->prev = temp;
    start->edge = h2;
    end->edge->prev = h2;
    h2->prev->next = h2;
    temp = end->edge;

    // Allocate the new face to each edge of the new polygon
    while (temp != h2)
    {
        temp->face = f2;
        temp = temp->next;
    }
}

// Check if the given vertex is a notch point for two adjacent edges.
// Returns true if the vertex is a notch point, otherwise false.
bool is_notch_point(Vertex *a, Vertex *b, Vertex *c)
{
    /** @brief This function checks if the given vertex is a notch point for two adjacent edges.
     * It computes the cross product of the vectors between the vertices a, b, and c.
     * If the cross product is positive, the vertex is a notch point.
     * @param a Pointer to the vertex a
     * @param b Pointer to the vertex b
     * @param c Pointer to the vertex c
     * @return bool True if the vertex is a notch point, otherwise false.
     */
    // Compute vectors between the vertices a, b, and c
    double x1 = b->x - a->x;
    double y1 = b->y - a->y;
    double x2 = c->x - b->x;
    double y2 = c->y - b->y;
    double cross_product = x1 * y2 - x2 * y1; // If the cross product is positive, b is a notch point
    if (cross_product > 0)
    {
        return true;
    }

    // Otherwise, b is not a notch point
    return false;
}
// Function to obtain the LPVS (Line of sight polygon visibility set) given a set of points L and a set of points P
void obtain_LPVS(vector<Vertex *> &L, vector<Vertex *> &P, vector<Vertex *> &LPVS)
{
    /** @brief This function obtains the LPVS (Line of sight polygon visibility set) given a set of points L and a set of points P.
     * It iterates over all points in P and checks if the current point is not in L and is a notch point.
     * If the current point is a notch point, it is added to the LPVS.
     * @param L Vector of pointers to vertices in L
     * @param P Vector of pointers to vertices in P
     * @param LPVS Vector of pointers to vertices in the LPVS
     * @return void
     */
    // Iterate over all points in P
    for (int i = 0; i < P.size(); i++)
    {
        // Check if the current point is not in L and is a notch point
        if (find(L.begin(), L.end(), P[i]) == L.end() && is_notch_point(P[(i - 1 + P.size()) % P.size()], P[i], P[(i + 1) % P.size()]))
        {
            // Add the point to the LPVS if it is a notch point
            LPVS.push_back(P[i]);
        }
    }
    return;
}

// Function to find the vertex with the maximum x-coordinate in a given list of vertices.
// Returns the maximum x-coordinate found.
double max_x(vector<Vertex *> &L)
{
    /** @brief This function finds the vertex with the maximum x-coordinate in a given list of vertices.
     * It iterates over all vertices in the list and updates the maximum x-coordinate found.
     * @param L Vector of pointers to vertices
     * @return double The maximum x-coordinate found.
     * @return void
     */
    double max = L[0]->x; // Assume the first vertex has the maximum x-coordinate.
    for (int i = 1; i < L.size(); i++)
    {
        if (L[i]->x > max) // If a vertex with a larger x-coordinate is found,
        {
            max = L[i]->x; // update the maximum x-coordinate.
        }
    }
    return max; // Return the maximum x-coordinate found.
}

// Find the vertex with the minimum x-coordinate in the given vector of vertices
double min_x(vector<Vertex *> &L)
{
    /** @brief This function finds the vertex with the minimum x-coordinate in a given list of vertices.
     * It iterates over all vertices in the list and updates the minimum x-coordinate found.
     * @param L Vector of pointers to vertices
     * @return double The minimum x-coordinate found.
     * @return void
     */
    // Set the minimum x-coordinate initially to be the x-coordinate of the first vertex in the vector
    double min = L[0]->x;
    // Iterate over all the remaining vertices in the vector
    for (int i = 1; i < L.size(); i++)
    {
        // If the x-coordinate of the current vertex is less than the current minimum x-coordinate, update the minimum
        if (L[i]->x < min)
        {
            min = L[i]->x;
        }
    }
    // Return the minimum x-coordinate
    return min;
}

// Find the maximum y-coordinate value among the vertices in the given list
double max_y(vector<Vertex *> &L)
{

    /** @brief This function finds the maximum y-coordinate value among the vertices in the given list.
     * It iterates over all vertices in the list and updates the maximum y-coordinate value found.
     * @param L Vector of pointers to vertices
     * @return double The maximum y-coordinate value found.
     */
    double max = L[0]->y;              // set the current maximum to the y-coordinate of the first vertex
    for (int i = 1; i < L.size(); i++) // iterate through the remaining vertices in the list
    {
        if (L[i]->y > max) // if the y-coordinate of the current vertex is greater than the current maximum
        {
            max = L[i]->y; // update the current maximum to the y-coordinate of the current vertex
        }
    }
    return max; // return the maximum y-coordinate value found
}

// Returns the minimum y coordinate of the vertices in the given vector
double min_y(vector<Vertex *> &L)
{
    /** @brief This function finds the minimum y-coordinate value among the vertices in the given list.
     * It iterates over all vertices in the list and updates the minimum y-coordinate value found.
     * @param L Vector of pointers to vertices
     * @return double The minimum y-coordinate value found.
     */
    // Set the minimum y-coordinate initially to be the y-coordinate of the first vertex in the vector
    double min = L[0]->y;
    // Iterate over all the remaining vertices in the vector
    for (int i = 1; i < L.size(); i++)
    {
        // If the y-coordinate of the current vertex is less than the current minimum y-coordinate, update the minimum
        if (L[i]->y < min)
        {
            min = L[i]->y;
        }
    }
    // Return the minimum y-coordinate
    return min;
}

// check vertex in rectangle
bool check_vertex_in_R(Vertex *v, double right, double left, double top, double bottom)
{
    /** @brief This function checks whether a vertex lies inside a rectangle or not.
     * @param v Pointer to the vertex to be checked
     * @param right The rightmost x-coordinate of the rectangle
     * @param left The leftmost x-coordinate of the rectangle
     * @param top The topmost y-coordinate of the rectangle
     * @param bottom The bottommost y-coordinate of the rectangle
     * @return bool True if the vertex lies inside the rectangle, false otherwise.
     */
    if (v->x >= left && v->x <= right && v->y >= bottom && v->y <= top)
    {
        return true;
    }
    return false;
}

// This function checks whether a point lies inside a polygon or not.
// s is the starting point of the polygon edge, n is a notch point on the edge,
// and p is the point to be checked.
bool check_inside_poly(Vertex *s, Vertex *n, Vertex *p)
{
    /** @brief This function checks whether a point lies inside a polygon or not.
     * @param s Pointer to the starting point of the polygon edge
     * @param n Pointer to a notch point on the edge
     * @param p Pointer to the point to be checked
     * @return bool True if the point lies inside the polygon, false otherwise.
     */
    // Calculate the vector from s to n (notch vector)
    double notch_vec_x = n->x - s->x; // a
    double notch_vec_y = n->y - s->y; // b

    // Calculate the vector from s to p (polygon vector)
    double poly_vec_x = p->x - s->x; // c
    double poly_vec_y = p->y - s->y; // d

    // Calculate the cross product of notch vector and polygon vector
    double cross_product = poly_vec_x * notch_vec_y - poly_vec_y * notch_vec_x;

    // If the cross product is greater than 0, the point lies inside the polygon
    if (cross_product > 0)
    {
        // Return true if the point is a notch point
        // (i.e., if it lies on an edge of the polygon)
        return true;
    }

    // Otherwise, the point lies outside the polygon
    return false;
}

// This function calculates the angle between three points: a, b, and c.
bool calculate_angle(Vertex *a, Vertex *b, Vertex *c)
{
    /** @brief This function calculates the angle between three points: a, b, and c.
     * @param a Pointer to the first point
     * @param b Pointer to the second point
     * @param c Pointer to the third point
     * @return bool True if the angle is reflex (greater than 180 degrees), false otherwise.
     */
    // Calculate the vectors from point a to b and from point b to c.
    double x1 = b->x - a->x;
    double y1 = b->y - a->y;
    double x2 = c->x - b->x;
    double y2 = c->y - b->y;
    // Calculate the cross product and dot product of the two vectors.
    // These are used to calculate the angle between the vectors.
    double cross_product = x1 * y2 - x2 * y1;
    double dot_product = x1 * x2 + y1 * y2;

    // Calculate the angle between the two vectors using the atan2 function.
    // This gives the angle in radians.
    double angle = atan2(cross_product, dot_product);

    // Convert the angle to degrees and check if it is greater than 180 degrees.
    // If it is, then the angle is reflex (greater than 180 degrees) and we return true.
    if (angle * 180 / 3.141592 > 180)
    {
        return true;
    }
    // Otherwise, the angle is not reflex and we return false.
    return false;
}

// This function determines if a given point is inside a polygon defined by its vertices.
int if_notch_in_poly(vector<Vertex *> vertices, double testx, double testy)
{
    /** @brief This function determines if a given point is inside a polygon defined by its vertices.
     * @param vertices Vector of pointers to vertices defining the polygon
     * @param testx The x-coordinate of the point to be checked
     * @param testy The y-coordinate of the point to be checked
     * @return int 1 if the point is inside the polygon, 0 if it is on the edge, and -1 if it is outside the polygon.
     */
    // Get the number of vertices in the polygon.
    int nvert = vertices.size();

    // Append the first vertex to the end of the vector and reverse the order of the vertices.
    vertices.push_back(vertices[0]);
    reverse(vertices.begin(), vertices.end());

    // Create separate vectors for the x-coordinates and y-coordinates of the vertices.
    vector<double> vertx;
    vector<double> verty;
    for (int i = 0; i < vertices.size(); i++)
    {
        vertx.push_back(vertices[i]->x);
        verty.push_back(vertices[i]->y);
    }

    // Initialize variables for the loop.
    int i, j, c = 0;
    for (i = 0, j = nvert - 1; i < nvert; j = i++)
    {
        // Determine if the test point's y-coordinate is within the current edge's y-coordinates.
        if (((verty[i] > testy) != (verty[j] > testy)) &&
            (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
            c = !c;
        // Increment or decrement the counter based on whether the current edge crosses the test point's x-axis.
    }

    // Return 0 if the point is outside the polygon, 1 if it is inside.
    return c;
}

// Implementing Algorithm 1
void solve(DCEL *polygon)
{
    /** @brief This function implements Algorithm 1 from the paper.
     * @param polygon Pointer to the DCEL representing the polygon
     * @return void
     */
    // Extract the vertices of the polygon into a vector P
    vector<Vertex *> P = polygon->vertices;

    // Create a new vector L and add the first vertex of P to it
    vector<Vertex *> L;
    L.push_back(P[0]);

    // While there are more than three vertices in P:
    while (P.size() > 3)
    {

        // Take the last vertex in L as v1
        Vertex *v1 = L[L.size() - 1];

        // Find the next vertex v2 in P that follows v1
        Vertex *v2;
        for (int i = 0; i < P.size(); i++)
        {
            if (P[i] == v1)
            {
                v2 = P[(i + 1) % P.size()];
                break;
            }
        }

        // Clear the contents of L and add v1 and v2 to it
        L.clear();
        L.push_back(v1);
        L.push_back(v2);

        // 3.3
        int ind;
        for (int i = 0; i < P.size(); i++)
        {

            if (P[i] == v2)
            {
                ind = i;
                break;
            }
        }
        ind = (ind + 1) % P.size();
        for (int i = ind; L.size() < P.size(); i = (i + 1) % P.size())
        {

            if (!is_notch_point(L[L.size() - 2], L[L.size() - 1], P[i]) &&
                !is_notch_point(L[L.size() - 1], P[i], L[0]) &&
                !is_notch_point(P[i], L[0], L[1]))
            {
                L.push_back(P[i]);
            }

            else
                break;
        }

        vector<Vertex *> LPVS;
        // 3.4 begins
        if (L.size() != P.size())
        {

            // 3.4.1 begins
            obtain_LPVS(L, P, LPVS);
            // 3.4.2 begins
            if (LPVS.size() > 0)
            {

                double right_edge = max_x(L);
                double left_edge = min_x(L);
                double top_edge = max_y(L);
                double bottom_edge = min_y(L);

                if (LPVS.size() > 0)
                {

                    // Repeat until
                    for (int i = 0; i < LPVS.size(); i++)
                    {

                        if (!check_vertex_in_R(LPVS[i], right_edge, left_edge, top_edge, bottom_edge))
                        {

                            LPVS.erase(LPVS.begin() + i);
                            i = -1;
                        }
                    }

                    // This is the VTR thing
                    for (int i = 0; i < LPVS.size(); i++)
                    {
                        Vertex *notch = LPVS[i];

                        // added this
                        if (!if_notch_in_poly(L, notch->x, notch->y))
                        {
                            continue;
                        }

                        if (!check_inside_poly(L[0], notch, L[L.size() - 1]))
                        {
                            continue;
                        }

                        for (int j = 1; j < L.size() && L.size() > 2; j++)
                        {
                            if (check_inside_poly(L[0], notch, L[j]))
                            {
                                L.erase(L.begin() + j);
                                j = 0;
                            }
                        }
                    }
                }
            }
        }
        // 3.5
        if (L.size() > 2 && L.size() != P.size())
        {

            // 3.5.1 begins
            split_face(L[0], L[L.size() - 1], polygon);

            // 3.5.2 begins
            for (int i = 1; i < L.size() - 1; i++)
            {
                for (int j = 0; j < P.size(); j++)
                {
                    if (L[i] == P[j])
                    {
                        P.erase(P.begin() + j);
                        j = -1;
                    }
                }
            }
        }
        if (L.size() == P.size())
        {
            P.clear();
        }
    }
}

// This function takes in three vertices p, q, and r as arguments.
// It computes the cross product of two vectors formed by the edges (p,q) and (q,r) using their 2D Cartesian coordinates.
// If the cross product is positive, it means the angle formed by these two edges is greater than 180 degrees, so it returns false.
// Otherwise, it returns true, indicating that the angle formed by these two edges is less than or equal to 180 degrees.
bool merge_angle(Vertex *p, Vertex *q, Vertex *r)
{
    /** @brief This function computes the cross product of two vectors formed by the edges (p,q) and (q,r) using their 2D Cartesian coordinates.
     * @param p The first vertex of the first edge.
     * @param q The second vertex of the first edge.
     * @param r The second vertex of the second edge.
     * @return Returns true if the angle formed by the two edges is less than or equal to 180 degrees, and false otherwise.
     */
    // Calculate the coordinates of the two vectors.
    double one_x = q->x - p->x;
    double one_y = q->y - p->y;
    double two_x = r->x - q->x;
    double two_y = r->y - q->y;

    // Compute the cross product of the two vectors.
    double prod = one_x * two_y - one_y * two_x;

    // If the cross product is positive, the angle is greater than 180 degrees.
    if (prod > 0)
    {
        return false;
    }

    // Otherwise, the angle is less than or equal to 180 degrees.
    return true;
}

// This function takes a DCEL (Doubly Connected Edge List) representing a polygon as an input,
// and merges any diagonals of the polygon that satisfy a certain condition.

void merge_polygon(DCEL *polygon)
{
    /** @brief This function merges any diagonals of the polygon that satisfy a certain condition.
     * @param polygon The DCEL representing the polygon.
     * @return Returns nothing.
     */
    // Iterate over all diagonals in the polygon.
    for (int i = 0; i < polygon->diagonals.size(); i++)
    {
        // Get the two HalfEdges associated with the diagonal.
        HalfEdge *h1 = polygon->diagonals[i];
        HalfEdge *h2 = h1->twin;

        // Get the Vertices associated with the diagonal.
        Vertex *b = h1->origin;
        Vertex *c = h1->prev->origin;
        Vertex *a = h2->next->next->origin;
        Vertex *d = h2->prev->origin;
        Vertex *e = h2->origin;
        Vertex *f = h1->next->next->origin;

        // Check if the diagonal satisfies the condition for merging.
        if (merge_angle(c, b, a) && merge_angle(d, e, f))
        {
            // If the diagonal satisfies the condition, remove the old faces by setting their is_valid flags to false.
            h1->face->is_valid = false;
            h2->face->is_valid = false;

            // Create a new face.
            Face *new_face = new Face();

            // Update the face pointers for all edges on both old faces to point to the new face.
            HalfEdge *temp = h1;
            while (temp->next != h1)
            {
                temp->face = new_face;
                temp = temp->next;
            }
            HalfEdge *temp2 = h2;
            while (temp2->next != h2)
            {
                temp2->face = new_face;
                temp2 = temp2->next;
            }

            // Add the new face to the list of faces in the polygon.
            polygon->faces.push_back(new_face);

            // Allocate an edge to the new face.
            new_face->edge = h1->next;

            // Perform edge updates to complete the merge.
            h1->prev->next = h2->next;
            h2->next->prev = h1->prev;
            h1->next->prev = h2->prev;
            h2->prev->next = h1->next;

            // Free the memory used by the old HalfEdges.
            delete (h1);
            delete (h2);
        }
        else
        {
            // If the diagonal does not satisfy the condition for merging, continue to the next diagonal.
            continue;
        }
    }
}

void merge(DCEL *polygon)
{
    int m = polygon->diagonals.size();
    int NP = m + 1;
    vector<HalfEdge *> LLe = polygon->diagonals;
    vector<bool> LDP(NP, true);
    vector<int> LUP(NP);
    map<Vertex *, vector<pair<int, Vertex *>>> LP;
    for (int i = 1; i <= NP; i++)
    {
        LUP[i] = i;
    }
    for (int j = 1; j <= m; j++)
    {
        Vertex *vs = LLe[j]->origin;
        Vertex *vt = LLe[j]->next->origin;
        if ((LP[vs].size() > 2 && LP[vt].size() > 2) || (LP[vs].size() > 2 && is_notch_point(vt->edge->prev->origin, vt, vt->edge->next->origin)) || (LP[vt].size() > 2 && is_notch_point(vs->edge->prev->origin, vs, vs->edge->next->origin)) || (is_notch_point(vt->edge->prev->origin, vt, vt->edge->next->origin) && is_notch_point(vs->edge->prev->origin, vs, vs->edge->next->origin)))
        {
            Vertex *j2 = vt;
            Vertex *i2 = vs;
            HalfEdge *j_edge = LLe[j];
            while (j_edge->origin != vt)
            {
                j_edge = j_edge->next;
            }
            Vertex *j3 = j_edge->next->origin;
            j_edge = LLe[j];
            while (j_edge->origin != vs)
            {
                j_edge = j_edge->next;
            }
            Vertex *i1 = j_edge->prev->origin;
            int u;
            for (auto it : LP[vt])
            {
                if (it.second == vs)
                {
                    u = it.first;
                    break;
                }
            }
            HalfEdge *u_edge = LLe[u];
            while (u_edge->origin != vt)
            {
                u_edge = u_edge->next;
            }
            Vertex *j1 = u_edge->prev->origin;
            u_edge = LLe[j];
            while (u_edge->origin != vs)
            {
                u_edge = u_edge->next;
            }
            Vertex *i3 = u_edge->next->origin;
            if (merge_angle(i1, i2, i3) && merge_angle(j1, j2, j3))
            {
                NP++;
                LDP.push_back(true);
                LDP[j] = false;
                LDP[u] = false;
                LUP[j] = NP;
                LUP[u] = NP;
                for (int h = 1; h < NP; h++)
                {
                    if (LUP[h] == j || LUP[h] == u)
                    {
                        LUP[h] = NP;
                    }
                }
            }
        }
    }
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    cout.tie(NULL);

    // write file pointers that will be used to read and write to files from the folder Scripts+Txt
    // freopen("../Scripts+TXT_Files/input.txt", "r", stdin);
    // freopen("../Scripts+TXT_Files/output.txt", "w+", stdout);

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w+", stdout);

    FILE *fp = fopen("merge.txt", "w+");
    int n;
    cin >> n;

    vector<pair<double, double>> points(n);
    for (int i = 0; i < n; i++)
    {
        cin >> points[i].first >> points[i].second;
    }
    //    reverse(points.begin(), points.end());
    DCEL *init_polygon = new DCEL();

    // Add clock
    clock_t start, end;
    start = clock();

    init_polygon->create_polygon(points);
    solve(init_polygon);

    cout << init_polygon->faces.size() << endl;
    for (int i = 0; i < init_polygon->faces.size(); i++)
    {
        if (init_polygon->faces[i]->is_valid)
            traverse_face(init_polygon->faces[i]);
    }

    merge_polygon(init_polygon);

    end = clock();

    // Print the time taken by the program
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(5);
    cout << " sec " << endl;

    int cunt = 0;
    for (int i = 0; i < init_polygon->faces.size(); i++)
    {
        if (init_polygon->faces[i]->is_valid)
            cunt++;
    }
    fprintf(fp, "%d\n", cunt);
    for (int i = 0; i < init_polygon->faces.size(); i++)
    {
        if (init_polygon->faces[i]->is_valid)
            traverse_face_merge(init_polygon->faces[i], fp);
    }
    return 0;
}