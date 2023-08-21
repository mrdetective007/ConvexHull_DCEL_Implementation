/** @file mainpage.h
 * @brief This file contains the analysis of the code.
 */
/** @mainpage Analysis of the code
 *
 * @authors Khushil Kataria 2020A7PS2086H
 * @authors Rahil Amish Sanghavi 2020A7PS2047H
 * @authors Jeet Shah 2020A7PS0155H
 * @authors Luv Ghilothia 2020A7PS1700H
 *
 * This is the main page for the implementation of Algorithm 1 from the paper "Algorithms for the
 * decomposition of a polygon into convex polygons" by J. Fernandez, L. Canovas, B. Pelegrin.
 * The algorithm takes as input a polygon represented as a set of points, and outputs a set of
 * convex polygons that cover the original polygon. The algorithm works by recursively splitting
 * the polygon into smaller convex pieces until the entire polygon is covered.
 * This implementation uses a divide-and-conquer approach to recursively split the polygon into
 * smaller sub-polygons. The algorithm first computes a diagonal of the polygon, which is a line
 * segment connecting two non-adjacent vertices that lies entirely inside the polygon. This diagonal
 * splits the polygon into two sub-polygons, which are then processed recursively. The algorithm
 * terminates when the sub-polygon is a triangle, which is already convex.
 * This implementation has a time complexity of O(n^3), where n is the number of vertices in the
 * polygon. The algorithm is guaranteed to terminate for any simple polygon with non-intersecting edges.
 * @section F Time Complexities
 * @subsection a void DCEL::create_polygon( vector<pair<double, double>> &arr )
 * O(n), where n is the number of vertices in the polygon.
 * @subsection b void traverse_face( Face *f )
 * O(n), where n is the number of vertices in the polygon.
 * @subsection c void traverse_face_merge( Face *f, FILE *fp )
 * O(n), where n is the number of vertices in the polygon.
 * @subsection d void split_face( Vertex *start, Vertex *end, DCEL *polygon )
 * O(n), where n is the number of vertices in the polygon.
 * @subsection e bool is_notch_point( Vertex *a, Vertex *b, Vertex *c )
 * O(1)
 * @subsection f void obtain_LPVS( vector<Vertex *> &L, vector<Vertex *> &P, vector<Vertex *> &LPVS )
 * O(nlogn), where n is the number of vertices in the polygon.
 * @subsection g double max_x( vector<Vertex *> &L )
 * O(1)
 * @subsection h double min_x( vector<Vertex *> &L )
 * O(1)
 * @subsection i double max_y( vector<Vertex *> &L )
 * O(1)
 * @subsection j double min_y( vector<Vertex *> &L )
 * O(1)
 * @subsection k bool check_vertex_in_R( Vertex *v, double right, double left, double top, double bottom )
 * O(1)
 * @subsection l bool check_inside_poly( Vertex *s, Vertex *n, Vertex *p )
 * O(1)
 * @subsection m int if_notch_in_poly( vector<Vertex *> vertices, double testx, double testy )
 * O(n), where n is the number of vertices in the polygon.
 * @subsection nm void solve( DCEL *polygon )
 * Ω(n²logn) and O(n³), where n is the number of vertices in the polygon.
 * @subsection o bool merge_angle( Vertex *p, Vertex *q, Vertex *r )
 * O(1)
 * @subsection p void merge_polygon( DCEL *polygon )
 * O(n²), where n is the number of vertices in the polygon.
 * @subsection q Overall Complexity
 * O(n³), where n is the number of vertices in the polygon.
 * @section Examples From the paper
 * @subsection Example_1
 * \image html images/example-1.png "Example 1" width=1000cm
 * @subsection Example_2
 * \image html images/example-2.png "Example 2" width=1000cm
 * @subsection Results
 * \image html images/result.png "Results" width=1000cm
 * @section Polygon_1
 * 5 points
 * Runtimes:
 * 1. 0.000206
 * 2. 0.000131
 * 3. 0.000109 <br>
 * <strong> Average runtime: 0.000149 </strong>
 * \image html images/input-5.png "Output" width=1000cm
 * @section Polygon_2
 * 9 points
 * Runtimes:
 * 1. 0.000159
 * 2. 0.000158
 * 3. 0.000169 <br>
 * <strong> Average runtime: 0.000162 </strong>
 * \image html images/input-9.png "Output" width=1000cm
 * @section Polygon_3
 * 10 points
 * Runtimes:
 * 1. 0.000246
 * 2. 0.000259
 * 3. 0.000266 <br>
 * <strong> Average runtime: 0.000257 </strong>
 * \image html images/input-10.png "Output" width=1000cm
 * @section Polygon_4
 * 13 points
 * Runtimes:
 * 1. 0.000205
 * 2. 0.000203
 * 3. 0.000208 <br>
 * <strong> Average runtime: 0.000205 </strong>
 * \image html images/input-13.png "Output" width=1000cm
 * @section Polygon_5
 * 15 points
 * Runtimes:
 * 1. 0.000245
 * 2. 0.000246
 * 3. 0.000242 <br>
 * <strong> Average runtime: 0.000244 </strong>
 * \image html images/input-15.png "Output" width=1000cm
 * @section Polygon_6
 * 20 points
 * Runtimes:
 * 1. 0.000346
 * 2. 0.000353
 * 3. 0.000357 <br>
 * <strong> Average runtime: 0.000352 </strong>
 * \image html images/input-20.png "Output" width=1000cm
 * @section Polygon_7
 * 22 points
 * Runtimes:
 * 1. 0.000280
 * 2. 0.000341
 * 3. 0.000334 <br>
 * <strong> Average runtime: 0.000318 </strong>
 * \image html images/input-22.png "Output" width=1000cm
 * @section Polygon_8
 * 25 points
 * Runtimes:
 * 1. 0.000394
 * 2. 0.000412
 * 3. 0.000400 <br>
 * <strong> Average runtime: 0.000402 </strong>
 * \image html images/input-25.png "Output" width=1000cm
 * @section Polygon_9
 * 37 points
 * Runtimes:
 * 1. 0.000684
 * 2. 0.000682
 * 3. 0.000694 <br>
 * <strong> Average runtime: 0.000687 </strong>
 * \image html images/input-37.png "Output" width=1000cm
 * @section Polygon_10
 * 39 points
 * Runtimes:
 * 1. 0.000680
 * 2. 0.000661
 * 3. 0.000689 <br>
 * <strong> Average runtime: 0.000676 </strong>
 * \image html images/input-39.png "Output" width=1000cm
 * @section Polygon_11
 * 50 points
 * 1. 0.001305
 * 2. 0.001310
 * 3. 0.001298 <br>
 * <strong> Average runtime: 0.001304 </strong>
 * \image html images/input-50.png "Output" width=1000cm
 * @section Polygon_12
 * 69 points
 * Runtimes:
 * 1. 0.001547
 * 2. 0.001520
 * 3. 0.001701 <br>
 * <strong> Average runtime: 0.001589 </strong>
 * \image html images/input-69.png "Output" width=1000cm
 * @section Polygon_13
 * 100 points
 * Runtimes:
 * 1. 0.004304
 * 2. 0.004322
 * 3. 0.004296 <br>
 * <strong> Average runtime: 0.004307 </strong>
 * \image html images/input-100.png "Output" width=1000cm
 * @section Polygon_14
 * 500 points
 * Runtimes:
 * 1. 0.058671
 * 2. 0.056839
 * 3. 0.057634 <br>
 * <strong> Average runtime: 0.057714 </strong>
 * \image html images/input-500.png "Output" width=1000cm
 * @section Polygon_15
 * 650 points
 * Runtimes:
 * 1. 0.094677
 * 2. 0.089513
 * 3. 0.089793 <br>
 * <strong> Average runtime: 0.091328 </strong>
 * \image html images/input-650.png "Output" width=1000cm
 * @section Polygon_16
 * 700 points
 * Runtimes:
 * 1. 0.111304
 * 2. 0.109460
 * 3. 0.109284 <br>
 * <strong> Average runtime: 0.110016 </strong>
 * \image html images/input-700.png "Output" width=1000cm
 * @section Polygon_17
 * 850 points
 * Runtimes:
 * 1. 0.163370
 * 2. 0.178663
 * 3. 0.165801
 * <strong> Average runtime: 0.169245 </strong>
 * \image html images/input-850.png "Output" width=1000cm
 * @section Polygon_18
 * 1000 points
 * Runtimes:
 * 1. 0.244714
 * 2. 0.240066
 * 3. 0.255304 <br>
 * <strong> Average runtime: 0.246695 </strong>
 * \image html images/input-1000.png "Output" width=1000cm
 * @section Polygon_19
 * 1280 points
 * Runtimes:
 * 1. 0.367157
 * 2. 0.361994
 * 3. 0.390300 <br>
 * <strong> Average runtime: 0.373150 </strong>
 * \image html images/input-1280.png "Output" width=1000cm
 * @section Polygon_20
 * 1410 points
 * Runtimes:
 * 1. 0.481227
 * 2. 0.487941
 * 3. 0.510249 <br>
 * <strong> Average runtime: 0.493139 </strong>
 * \image html images/input-1410.png "Output" width=1000cm
 * @section Polygon_21
 * 1665 points
 * Runtimes:
 * 1. 0.578529
 * 2. 0.577223
 * 3. 0.624089 <br>
 * <strong> Average runtime: 0.593280 </strong>
 * \image html images/input-1665.png "Output" width=1000cm
 * @section Polygon_22
 * 1873 points
 * Runtimes:
 * 1. 1.042004
 * 2. 0.929109
 * 3. 0.912022 <br>
 * <strong> Average runtime: 0.961045 </strong>
 * \image html images/input-1873.png "Output" width=1000cm
 * @section Polygon_23
 * 2000 points
 * Runtimes:
 * 1. 0.881860
 * 2. 0.877704
 * 3. 0.876115 <br>
 * <strong> Average runtime: 0.878560 </strong>
 * \image html images/input-2000.png "Output" width=1000cm
 * @section Polygon_24
 * 2650 points
 * Runtimes:
 * 1. 1.772407
 * 2. 1.777430
 * 3. 1.771009
 * <strong> Average runtime: 1.773615 </strong>
 * \image html images/input-2650.png "Output" width=1000cm
 * @section Polygon_25
 * 3200 points
 * Runtimes:
 * 1. 2.664388
 * 2. 2.662064
 * 3. 2.667416 <br>
 * <strong> Average runtime: 2.664623 </strong>
 * \image html images/input-3200.png "Output" width=1000cm
 * @section Polygon_26
 * 3800 points
 * Runtimes:
 * 1. 2.697847
 * 2. 2.660535
 * 3. 2.673392 <br>
 * <strong> Average runtime: 2.677258 </strong>
 * \image html images/input-3800.png "Output" width=1000cm
 * @section Polygon_27
 * 4267 points
 * Runtimes:
 * 1. 4.975617
 * 2. 4.857599
 * 3. 4.829600 <br>
 * <strong> Average runtime: 4.887605 </strong>
 * \image html images/input-4267.png "Output" width=1000cm
 * @section Polygon_28
 * 5000 points
 * Runtimes:
 * 1. 6.646950
 * 2. 6.633213
 * 3. 6.701837 <br>
 * <strong> Average runtime: 6.660667</strong>
 * \image html images/input-5000.png "Output" width=1000cm
 */