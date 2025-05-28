#ifndef FIELDS_H
#define FIELDS_H

#include "mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Data structures for field variables */

typedef struct {
    double ux, uy;  /* Displacement components */
} Displacement;

typedef struct {
    double xx, yy, xy;  /* Stress tensor components (2D) */
} Stress;

typedef struct {
    double xx, yy, xy;  /* Strain tensor components (2D) */
} Strain;

typedef struct {
    int node_id;
    Point2D coord;
    Displacement disp;
    Stress stress;
    Strain strain;
    double strain_energy_density;  /* W = 0.5 * sigma_ij * epsilon_ij */
} FieldData;

typedef struct {
    int num_points;
    FieldData *data;
    Mesh *mesh;
} FieldSet;

/* Function prototypes */

/**
 * Create a field set for a given mesh
 * @param mesh Pointer to mesh
 * @return Pointer to created field set
 */
FieldSet* create_field_set(Mesh *mesh);

/**
 * Generate synthetic stress field (polynomial approximation)
 * Useful for testing and validation
 * @param fields Pointer to field set
 * @param crack_tip Crack tip location
 * @param applied_stress Applied far-field stress
 */
void generate_synthetic_stress_field(FieldSet *fields, Point2D crack_tip, 
                                   double applied_stress);

/**
 * Generate Mode I crack field using Williams series expansion
 * @param fields Pointer to field set
 * @param crack_tip Crack tip location
 * @param K_I Mode I stress intensity factor
 * @param num_terms Number of terms in Williams series
 */
void generate_williams_field(FieldSet *fields, Point2D crack_tip, 
                           double K_I, int num_terms);

/**
 * Load field data from file
 * @param filename Input file name
 * @param mesh Associated mesh
 * @return Pointer to loaded field set
 */
FieldSet* load_fields_from_file(const char *filename, Mesh *mesh);

/**
 * Save field data to file
 * @param fields Pointer to field set
 * @param filename Output file name
 */
void save_fields_to_file(const FieldSet *fields, const char *filename);

/**
 * Interpolate field values at a given point using bilinear interpolation
 * @param fields Pointer to field set
 * @param point Point where to interpolate
 * @param disp Interpolated displacement (output)
 * @param stress Interpolated stress (output)
 * @param strain_energy_density Interpolated strain energy density (output)
 * @return 0 on success, -1 on failure
 */
int interpolate_fields(const FieldSet *fields, Point2D point,
                      Displacement *disp, Stress *stress, 
                      double *strain_energy_density);

/**
 * Calculate strain from displacement gradients
 * @param du_dx Displacement gradient du/dx
 * @param du_dy Displacement gradient du/dy
 * @param dv_dx Displacement gradient dv/dx
 * @param dv_dy Displacement gradient dv/dy
 * @return Strain tensor
 */
Strain calculate_strain_from_gradients(double du_dx, double du_dy,
                                     double dv_dx, double dv_dy);

/**
 * Calculate strain energy density from stress and strain
 * @param stress Stress tensor
 * @param strain Strain tensor
 * @return Strain energy density W
 */
double calculate_strain_energy_density(const Stress *stress, const Strain *strain);

/**
 * Calculate displacement gradients at a point
 * @param fields Pointer to field set
 * @param point Point where to calculate gradients
 * @param du_dx Gradient du/dx (output)
 * @param du_dy Gradient du/dy (output)
 * @param dv_dx Gradient dv/dx (output)
 * @param dv_dy Gradient dv/dy (output)
 * @return 0 on success, -1 on failure
 */
int calculate_displacement_gradients(const FieldSet *fields, Point2D point,
                                   double *du_dx, double *du_dy,
                                   double *dv_dx, double *dv_dy);

/**
 * Update strain energy density for all field points
 * @param fields Pointer to field set
 */
void update_strain_energy_density(FieldSet *fields);

/**
 * Free field set memory
 * @param fields Pointer to field set to free
 */
void free_field_set(FieldSet *fields);

/**
 * Print field information at a specific point
 * @param fields Pointer to field set
 * @param point_id Point index
 */
void print_field_info(const FieldSet *fields, int point_id);

/**
 * Validate field data consistency
 * @param fields Pointer to field set
 * @return 0 if valid, error code otherwise
 */
int validate_field_data(const FieldSet *fields);

#endif /* FIELDS_H */ 