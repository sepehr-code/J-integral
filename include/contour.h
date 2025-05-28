#ifndef CONTOUR_H
#define CONTOUR_H

#include "mesh.h"
#include "fields.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Data structures for contour integration */

typedef struct {
    Point2D point;      /* Integration point coordinates */
    Point2D normal;     /* Outward normal vector */
    double ds;          /* Differential arc length */
    double weight;      /* Integration weight */
} ContourPoint;

typedef struct {
    int num_points;           /* Number of integration points */
    ContourPoint *points;     /* Array of integration points */
    Point2D center;           /* Contour center (usually crack tip) */
    double radius;            /* Contour radius */
    int contour_type;         /* 0: circular, 1: rectangular, 2: custom */
    char name[64];            /* Contour identifier */
} Contour;

typedef struct {
    int num_contours;         /* Number of contours */
    Contour *contours;        /* Array of contours */
    Point2D crack_tip;        /* Crack tip location */
} ContourSet;

/* Contour types */
#define CONTOUR_CIRCULAR    0
#define CONTOUR_RECTANGULAR 1
#define CONTOUR_CUSTOM      2

/* Function prototypes */

/**
 * Create a contour set around a crack tip
 * @param crack_tip Crack tip location
 * @param num_contours Number of contours to create
 * @return Pointer to created contour set
 */
ContourSet* create_contour_set(Point2D crack_tip, int num_contours);

/**
 * Add a circular contour around crack tip
 * @param contour_set Pointer to contour set
 * @param radius Contour radius
 * @param num_points Number of integration points
 * @param name Contour name
 * @return 0 on success, -1 on failure
 */
int add_circular_contour(ContourSet *contour_set, double radius, 
                        int num_points, const char *name);

/**
 * Add a rectangular contour around crack tip
 * @param contour_set Pointer to contour set
 * @param width Rectangle width
 * @param height Rectangle height
 * @param num_points_per_side Number of points per side
 * @param name Contour name
 * @return 0 on success, -1 on failure
 */
int add_rectangular_contour(ContourSet *contour_set, double width, double height,
                           int num_points_per_side, const char *name);

/**
 * Add a custom contour from point array
 * @param contour_set Pointer to contour set
 * @param points Array of contour points
 * @param num_points Number of points
 * @param name Contour name
 * @return 0 on success, -1 on failure
 */
int add_custom_contour(ContourSet *contour_set, Point2D *points, 
                      int num_points, const char *name);

/**
 * Calculate outward normal vectors for contour points
 * @param contour Pointer to contour
 */
void calculate_contour_normals(Contour *contour);

/**
 * Calculate differential arc lengths for contour segments
 * @param contour Pointer to contour
 */
void calculate_contour_arc_lengths(Contour *contour);

/**
 * Verify that contour is counter-clockwise oriented
 * @param contour Pointer to contour
 * @return 1 if counter-clockwise, 0 if clockwise, -1 on error
 */
int verify_contour_orientation(const Contour *contour);

/**
 * Reverse contour orientation
 * @param contour Pointer to contour
 */
void reverse_contour_orientation(Contour *contour);

/**
 * Check if a point is inside a contour
 * @param contour Pointer to contour
 * @param point Point to check
 * @return 1 if inside, 0 if outside
 */
int point_inside_contour(const Contour *contour, Point2D point);

/**
 * Find the closest contour point to a given location
 * @param contour Pointer to contour
 * @param point Reference point
 * @return Index of closest contour point
 */
int find_closest_contour_point(const Contour *contour, Point2D point);

/**
 * Refine contour by adding more integration points
 * @param contour Pointer to contour
 * @param refinement_factor Factor by which to increase points
 * @return 0 on success, -1 on failure
 */
int refine_contour(Contour *contour, int refinement_factor);

/**
 * Load contour set from file
 * @param filename Input file name
 * @param crack_tip Crack tip location
 * @return Pointer to loaded contour set
 */
ContourSet* load_contour_set_from_file(const char *filename, Point2D crack_tip);

/**
 * Save contour set to file
 * @param contour_set Pointer to contour set
 * @param filename Output file name
 */
void save_contour_set_to_file(const ContourSet *contour_set, const char *filename);

/**
 * Print contour information
 * @param contour Pointer to contour
 */
void print_contour_info(const Contour *contour);

/**
 * Print contour set information
 * @param contour_set Pointer to contour set
 */
void print_contour_set_info(const ContourSet *contour_set);

/**
 * Free contour memory
 * @param contour Pointer to contour to free
 */
void free_contour(Contour *contour);

/**
 * Free contour set memory
 * @param contour_set Pointer to contour set to free
 */
void free_contour_set(ContourSet *contour_set);

/**
 * Validate contour geometry
 * @param contour Pointer to contour
 * @return 0 if valid, error code otherwise
 */
int validate_contour(const Contour *contour);

#endif /* CONTOUR_H */ 