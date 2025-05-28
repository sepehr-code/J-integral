#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Data structures for 2D mesh and crack geometry */

typedef struct {
    double x, y;
} Point2D;

typedef struct {
    int id;
    Point2D coord;
} Node;

typedef struct {
    int nodes[4];  /* Quadrilateral element with 4 nodes */
    double area;
} Element;

typedef struct {
    Point2D start;
    Point2D end;
    Point2D tip;     /* Crack tip location */
    double length;
    double angle;    /* Crack orientation angle */
} Crack;

typedef struct {
    int nx, ny;           /* Number of nodes in x and y directions */
    int num_nodes;        /* Total number of nodes */
    int num_elements;     /* Total number of elements */
    Node *nodes;          /* Array of nodes */
    Element *elements;    /* Array of elements */
    Crack crack;          /* Crack geometry */
    double width, height; /* Domain dimensions */
} Mesh;

/* Function prototypes */

/**
 * Create a rectangular mesh with specified dimensions
 * @param nx Number of nodes in x-direction
 * @param ny Number of nodes in y-direction
 * @param width Domain width
 * @param height Domain height
 * @return Pointer to created mesh
 */
Mesh* create_rectangular_mesh(int nx, int ny, double width, double height);

/**
 * Define a crack in the mesh
 * @param mesh Pointer to mesh
 * @param start_x Crack start x-coordinate
 * @param start_y Crack start y-coordinate
 * @param end_x Crack end x-coordinate
 * @param end_y Crack end y-coordinate
 */
void define_crack(Mesh *mesh, double start_x, double start_y, 
                  double end_x, double end_y);

/**
 * Load mesh from file
 * @param filename Input file name
 * @return Pointer to loaded mesh
 */
Mesh* load_mesh_from_file(const char *filename);

/**
 * Save mesh to file
 * @param mesh Pointer to mesh
 * @param filename Output file name
 */
void save_mesh_to_file(const Mesh *mesh, const char *filename);

/**
 * Find the element containing a given point
 * @param mesh Pointer to mesh
 * @param point Point to locate
 * @return Element index (-1 if not found)
 */
int find_element_containing_point(const Mesh *mesh, Point2D point);

/**
 * Get node coordinates by index
 * @param mesh Pointer to mesh
 * @param node_id Node index
 * @return Node coordinates
 */
Point2D get_node_coordinates(const Mesh *mesh, int node_id);

/**
 * Calculate distance between two points
 * @param p1 First point
 * @param p2 Second point
 * @return Distance
 */
double point_distance(Point2D p1, Point2D p2);

/**
 * Free mesh memory
 * @param mesh Pointer to mesh to free
 */
void free_mesh(Mesh *mesh);

/**
 * Print mesh information
 * @param mesh Pointer to mesh
 */
void print_mesh_info(const Mesh *mesh);

#endif /* MESH_H */ 