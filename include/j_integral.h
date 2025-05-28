#ifndef J_INTEGRAL_H
#define J_INTEGRAL_H

#include "mesh.h"
#include "fields.h"
#include "contour.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Data structures for J-integral calculation */

typedef struct {
    double j_value;           /* J-integral value */
    double energy_term;       /* W * delta_1j * n_j term */
    double stress_term;       /* sigma_ij * du_i/dx_1 * n_j term */
    double contour_length;    /* Total contour length */
    int num_integration_points; /* Number of integration points used */
    char contour_name[64];    /* Associated contour name */
    double computation_time;  /* Time taken for calculation (seconds) */
} JIntegralResult;

typedef struct {
    int num_results;          /* Number of J-integral results */
    JIntegralResult *results; /* Array of results for different contours */
    Point2D crack_tip;        /* Crack tip location */
    double mean_j_value;      /* Mean J-integral value */
    double std_deviation;     /* Standard deviation of J values */
    double path_independence_error; /* Measure of path independence */
    char analysis_name[128];  /* Analysis identifier */
} JIntegralAnalysis;

/* Integration methods */
#define INTEGRATION_MIDPOINT    0
#define INTEGRATION_TRAPEZOIDAL 1
#define INTEGRATION_SIMPSON     2
#define INTEGRATION_GAUSS       3

/* Function prototypes */

/**
 * Calculate J-integral for a single contour
 * @param fields Pointer to field set
 * @param contour Pointer to contour
 * @param integration_method Integration method to use
 * @return J-integral result structure
 */
JIntegralResult calculate_j_integral_single_contour(const FieldSet *fields,
                                                   const Contour *contour,
                                                   int integration_method);

/**
 * Calculate J-integral for multiple contours (path independence study)
 * @param fields Pointer to field set
 * @param contour_set Pointer to contour set
 * @param integration_method Integration method to use
 * @param analysis_name Name for this analysis
 * @return J-integral analysis structure
 */
JIntegralAnalysis* calculate_j_integral_multiple_contours(const FieldSet *fields,
                                                         const ContourSet *contour_set,
                                                         int integration_method,
                                                         const char *analysis_name);

/**
 * Calculate J-integral using midpoint rule
 * @param fields Pointer to field set
 * @param contour Pointer to contour
 * @return J-integral value
 */
double j_integral_midpoint_rule(const FieldSet *fields, const Contour *contour);

/**
 * Calculate J-integral using trapezoidal rule
 * @param fields Pointer to field set
 * @param contour Pointer to contour
 * @return J-integral value
 */
double j_integral_trapezoidal_rule(const FieldSet *fields, const Contour *contour);

/**
 * Calculate J-integral using Simpson's rule
 * @param fields Pointer to field set
 * @param contour Pointer to contour
 * @return J-integral value
 */
double j_integral_simpson_rule(const FieldSet *fields, const Contour *contour);

/**
 * Calculate the integrand at a specific contour point
 * J = ∫(W*δ₁ⱼ - σᵢⱼ*∂uᵢ/∂x₁)*nⱼ ds
 * @param fields Pointer to field set
 * @param point Contour point
 * @param integrand_value Calculated integrand value (output)
 * @param energy_term Energy term contribution (output)
 * @param stress_term Stress term contribution (output)
 * @return 0 on success, -1 on failure
 */
int calculate_integrand_at_point(const FieldSet *fields, 
                                const ContourPoint *point,
                                double *integrand_value,
                                double *energy_term,
                                double *stress_term);

/**
 * Perform convergence study with mesh refinement
 * @param base_mesh Base mesh for refinement
 * @param base_fields Base field set
 * @param crack_tip Crack tip location
 * @param max_refinement_levels Maximum refinement levels
 * @param contour_radius Contour radius for study
 * @return Array of J-integral values for different refinement levels
 */
double* convergence_study_mesh_refinement(const Mesh *base_mesh,
                                         const FieldSet *base_fields,
                                         Point2D crack_tip,
                                         int max_refinement_levels,
                                         double contour_radius);

/**
 * Perform convergence study with contour refinement
 * @param fields Pointer to field set
 * @param contour Base contour for refinement
 * @param max_refinement_levels Maximum refinement levels
 * @return Array of J-integral values for different contour refinements
 */
double* convergence_study_contour_refinement(const FieldSet *fields,
                                           const Contour *contour,
                                           int max_refinement_levels);

/**
 * Calculate path independence error
 * @param analysis Pointer to J-integral analysis
 * @return Path independence error measure
 */
double calculate_path_independence_error(const JIntegralAnalysis *analysis);

/**
 * Validate J-integral calculation
 * @param analysis Pointer to J-integral analysis
 * @param tolerance Acceptable tolerance for path independence
 * @return 0 if valid, error code otherwise
 */
int validate_j_integral_calculation(const JIntegralAnalysis *analysis, 
                                   double tolerance);

/**
 * Export J-integral results to file
 * @param analysis Pointer to J-integral analysis
 * @param filename Output file name
 * @param format Output format (0: text, 1: CSV, 2: binary)
 */
void export_j_integral_results(const JIntegralAnalysis *analysis,
                              const char *filename, int format);

/**
 * Load J-integral results from file
 * @param filename Input file name
 * @return Pointer to loaded J-integral analysis
 */
JIntegralAnalysis* load_j_integral_results(const char *filename);

/**
 * Print J-integral result summary
 * @param result Pointer to J-integral result
 */
void print_j_integral_result(const JIntegralResult *result);

/**
 * Print J-integral analysis summary
 * @param analysis Pointer to J-integral analysis
 */
void print_j_integral_analysis(const JIntegralAnalysis *analysis);

/**
 * Calculate statistics for J-integral analysis
 * @param analysis Pointer to J-integral analysis
 */
void calculate_j_integral_statistics(JIntegralAnalysis *analysis);

/**
 * Free J-integral result memory
 * @param result Pointer to result to free
 */
void free_j_integral_result(JIntegralResult *result);

/**
 * Free J-integral analysis memory
 * @param analysis Pointer to analysis to free
 */
void free_j_integral_analysis(JIntegralAnalysis *analysis);

/**
 * Compare J-integral with analytical solution (if available)
 * @param calculated_j Calculated J-integral value
 * @param analytical_j Analytical J-integral value
 * @param relative_error Relative error (output)
 * @return 0 on success, -1 on failure
 */
int compare_with_analytical_solution(double calculated_j, double analytical_j,
                                    double *relative_error);

/**
 * Calculate J-integral for Mode I crack using analytical formula
 * For validation purposes: J = K_I^2 / E' where E' = E for plane stress, E/(1-ν²) for plane strain
 * @param K_I Mode I stress intensity factor
 * @param E Young's modulus
 * @param nu Poisson's ratio
 * @param plane_stress 1 for plane stress, 0 for plane strain
 * @return Analytical J-integral value
 */
double analytical_j_integral_mode_I(double K_I, double E, double nu, int plane_stress);

#endif /* J_INTEGRAL_H */ 