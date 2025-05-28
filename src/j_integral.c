#include "../include/j_integral.h"
#include <string.h>
#include <time.h>

JIntegralResult calculate_j_integral_single_contour(const FieldSet *fields,
                                                   const Contour *contour,
                                                   int integration_method) {
    JIntegralResult result;
    memset(&result, 0, sizeof(JIntegralResult));
    
    if (!fields || !contour) {
        fprintf(stderr, "Error: Null pointer in J-integral calculation\n");
        return result;
    }
    
    clock_t start_time = clock();
    
    // Copy contour name
    strncpy(result.contour_name, contour->name, sizeof(result.contour_name) - 1);
    result.contour_name[sizeof(result.contour_name) - 1] = '\0';
    
    result.num_integration_points = contour->num_points;
    
    // Calculate total contour length
    result.contour_length = 0.0;
    for (int i = 0; i < contour->num_points; i++) {
        result.contour_length += contour->points[i].ds;
    }
    
    // Choose integration method
    switch (integration_method) {
        case INTEGRATION_MIDPOINT:
            result.j_value = j_integral_midpoint_rule(fields, contour);
            break;
        case INTEGRATION_TRAPEZOIDAL:
            result.j_value = j_integral_trapezoidal_rule(fields, contour);
            break;
        case INTEGRATION_SIMPSON:
            result.j_value = j_integral_simpson_rule(fields, contour);
            break;
        default:
            result.j_value = j_integral_midpoint_rule(fields, contour);
            break;
    }
    
    // Calculate individual terms for analysis
    result.energy_term = 0.0;
    result.stress_term = 0.0;
    
    for (int i = 0; i < contour->num_points; i++) {
        double integrand, energy_term, stress_term;
        if (calculate_integrand_at_point(fields, &contour->points[i], 
                                       &integrand, &energy_term, &stress_term) == 0) {
            result.energy_term += energy_term * contour->points[i].ds;
            result.stress_term += stress_term * contour->points[i].ds;
        }
    }
    
    clock_t end_time = clock();
    result.computation_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    return result;
}

JIntegralAnalysis* calculate_j_integral_multiple_contours(const FieldSet *fields,
                                                         const ContourSet *contour_set,
                                                         int integration_method,
                                                         const char *analysis_name) {
    if (!fields || !contour_set || contour_set->num_contours == 0) {
        fprintf(stderr, "Error: Invalid parameters for multiple contour analysis\n");
        return NULL;
    }
    
    JIntegralAnalysis *analysis = (JIntegralAnalysis*)malloc(sizeof(JIntegralAnalysis));
    if (!analysis) {
        fprintf(stderr, "Error: Memory allocation failed for J-integral analysis\n");
        return NULL;
    }
    
    analysis->num_results = contour_set->num_contours;
    analysis->crack_tip = contour_set->crack_tip;
    strncpy(analysis->analysis_name, analysis_name, sizeof(analysis->analysis_name) - 1);
    analysis->analysis_name[sizeof(analysis->analysis_name) - 1] = '\0';
    
    analysis->results = (JIntegralResult*)malloc(analysis->num_results * sizeof(JIntegralResult));
    if (!analysis->results) {
        fprintf(stderr, "Error: Memory allocation failed for J-integral results\n");
        free(analysis);
        return NULL;
    }
    
    // Calculate J-integral for each contour
    for (int i = 0; i < contour_set->num_contours; i++) {
        analysis->results[i] = calculate_j_integral_single_contour(fields, 
                                                                  &contour_set->contours[i],
                                                                  integration_method);
    }
    
    // Calculate statistics
    calculate_j_integral_statistics(analysis);
    
    return analysis;
}

double j_integral_midpoint_rule(const FieldSet *fields, const Contour *contour) {
    if (!fields || !contour) return 0.0;
    
    double j_value = 0.0;
    
    for (int i = 0; i < contour->num_points; i++) {
        double integrand, energy_term, stress_term;
        
        if (calculate_integrand_at_point(fields, &contour->points[i], 
                                       &integrand, &energy_term, &stress_term) == 0) {
            j_value += integrand * contour->points[i].ds;
        }
    }
    
    return j_value;
}

double j_integral_trapezoidal_rule(const FieldSet *fields, const Contour *contour) {
    if (!fields || !contour) return 0.0;
    
    double j_value = 0.0;
    
    for (int i = 0; i < contour->num_points; i++) {
        int next = (i + 1) % contour->num_points;
        
        double integrand1, integrand2;
        double energy_term, stress_term;
        
        // Calculate integrand at current and next points
        if (calculate_integrand_at_point(fields, &contour->points[i], 
                                       &integrand1, &energy_term, &stress_term) == 0 &&
            calculate_integrand_at_point(fields, &contour->points[next], 
                                       &integrand2, &energy_term, &stress_term) == 0) {
            
            // Trapezoidal rule: (f1 + f2) * ds / 2
            j_value += 0.5 * (integrand1 + integrand2) * contour->points[i].ds;
        }
    }
    
    return j_value;
}

double j_integral_simpson_rule(const FieldSet *fields, const Contour *contour) {
    if (!fields || !contour || contour->num_points < 3) return 0.0;
    
    double j_value = 0.0;
    
    // Simpson's rule requires even number of intervals
    int n = contour->num_points;
    if (n % 2 == 0) n--; // Use n-1 points if even
    
    for (int i = 0; i < n - 2; i += 2) {
        double integrand1, integrand2, integrand3;
        double energy_term, stress_term;
        
        // Calculate integrand at three consecutive points
        if (calculate_integrand_at_point(fields, &contour->points[i], 
                                       &integrand1, &energy_term, &stress_term) == 0 &&
            calculate_integrand_at_point(fields, &contour->points[i + 1], 
                                       &integrand2, &energy_term, &stress_term) == 0 &&
            calculate_integrand_at_point(fields, &contour->points[i + 2], 
                                       &integrand3, &energy_term, &stress_term) == 0) {
            
            // Simpson's rule: (f1 + 4*f2 + f3) * h / 3
            double h = (contour->points[i].ds + contour->points[i + 1].ds) / 2.0;
            j_value += (integrand1 + 4.0 * integrand2 + integrand3) * h / 3.0;
        }
    }
    
    return j_value;
}

int calculate_integrand_at_point(const FieldSet *fields, 
                                const ContourPoint *point,
                                double *integrand_value,
                                double *energy_term,
                                double *stress_term) {
    if (!fields || !point || !integrand_value || !energy_term || !stress_term) {
        return -1;
    }
    
    // Interpolate field values at the contour point
    Displacement disp;
    Stress stress;
    double strain_energy_density;
    
    if (interpolate_fields(fields, point->point, &disp, &stress, &strain_energy_density) != 0) {
        return -1;
    }
    
    // Calculate displacement gradients
    double du_dx, du_dy, dv_dx, dv_dy;
    if (calculate_displacement_gradients(fields, point->point, &du_dx, &du_dy, &dv_dx, &dv_dy) != 0) {
        return -1;
    }
    
    // J-integral formula: J = ∫(W*δ₁ⱼ - σᵢⱼ*∂uᵢ/∂x₁)*nⱼ ds
    // Where x₁ is the crack propagation direction (x-direction)
    
    // Energy term: W * δ₁ⱼ * nⱼ = W * n₁ (only j=1 contributes due to Kronecker delta)
    *energy_term = strain_energy_density * point->normal.x;
    
    // Stress term: -σᵢⱼ * ∂uᵢ/∂x₁ * nⱼ
    // = -(σ₁₁*∂u₁/∂x₁*n₁ + σ₁₂*∂u₁/∂x₁*n₂ + σ₂₁*∂u₂/∂x₁*n₁ + σ₂₂*∂u₂/∂x₁*n₂)
    // = -(σₓₓ*∂uₓ/∂x*nₓ + σₓᵧ*∂uₓ/∂x*nᵧ + σₓᵧ*∂uᵧ/∂x*nₓ + σᵧᵧ*∂uᵧ/∂x*nᵧ)
    
    *stress_term = -(stress.xx * du_dx * point->normal.x +
                     stress.xy * du_dx * point->normal.y +
                     stress.xy * dv_dx * point->normal.x +
                     stress.yy * dv_dx * point->normal.y);
    
    // Total integrand
    *integrand_value = *energy_term + *stress_term;
    
    return 0;
}

double* convergence_study_mesh_refinement(const Mesh *base_mesh,
                                         const FieldSet *base_fields,
                                         Point2D crack_tip,
                                         int max_refinement_levels,
                                         double contour_radius) {
    if (!base_mesh || !base_fields || max_refinement_levels <= 0) {
        return NULL;
    }
    
    double *j_values = (double*)malloc(max_refinement_levels * sizeof(double));
    if (!j_values) {
        fprintf(stderr, "Error: Memory allocation failed for convergence study\n");
        return NULL;
    }
    
    printf("Starting mesh refinement convergence study...\n");
    
    for (int level = 0; level < max_refinement_levels; level++) {
        int refinement_factor = 1 << level; // 1, 2, 4, 8, ...
        int nx_refined = (base_mesh->nx - 1) * refinement_factor + 1;
        int ny_refined = (base_mesh->ny - 1) * refinement_factor + 1;
        
        // Create refined mesh
        Mesh *refined_mesh = create_rectangular_mesh(nx_refined, ny_refined,
                                                    base_mesh->width, base_mesh->height);
        if (!refined_mesh) {
            j_values[level] = 0.0;
            continue;
        }
        
        // Copy crack definition
        refined_mesh->crack = base_mesh->crack;
        
        // Create field set for refined mesh (using synthetic fields for simplicity)
        FieldSet *refined_fields = create_field_set(refined_mesh);
        if (!refined_fields) {
            free_mesh(refined_mesh);
            j_values[level] = 0.0;
            continue;
        }
        
        // Generate synthetic stress field
        generate_synthetic_stress_field(refined_fields, crack_tip, 100.0);
        
        // Create contour set
        ContourSet *contour_set = create_contour_set(crack_tip, 1);
        add_circular_contour(contour_set, contour_radius, 64, "convergence_contour");
        
        // Calculate J-integral
        JIntegralResult result = calculate_j_integral_single_contour(refined_fields,
                                                                   &contour_set->contours[0],
                                                                   INTEGRATION_MIDPOINT);
        j_values[level] = result.j_value;
        
        printf("Level %d: %dx%d mesh, J = %.6e\n", level, nx_refined, ny_refined, j_values[level]);
        
        // Cleanup
        free_field_set(refined_fields);
        free_mesh(refined_mesh);
        free_contour_set(contour_set);
    }
    
    return j_values;
}

double* convergence_study_contour_refinement(const FieldSet *fields,
                                           const Contour *contour,
                                           int max_refinement_levels) {
    if (!fields || !contour || max_refinement_levels <= 0) {
        return NULL;
    }
    
    double *j_values = (double*)malloc(max_refinement_levels * sizeof(double));
    if (!j_values) {
        fprintf(stderr, "Error: Memory allocation failed for contour refinement study\n");
        return NULL;
    }
    
    printf("Starting contour refinement convergence study...\n");
    
    // Create a copy of the contour for refinement
    Contour refined_contour = *contour;
    refined_contour.points = (ContourPoint*)malloc(contour->num_points * sizeof(ContourPoint));
    memcpy(refined_contour.points, contour->points, contour->num_points * sizeof(ContourPoint));
    
    for (int level = 0; level < max_refinement_levels; level++) {
        // Calculate J-integral with current refinement
        JIntegralResult result = calculate_j_integral_single_contour(fields, &refined_contour,
                                                                   INTEGRATION_MIDPOINT);
        j_values[level] = result.j_value;
        
        printf("Level %d: %d points, J = %.6e\n", level, refined_contour.num_points, j_values[level]);
        
        // Refine contour for next iteration (except for last iteration)
        if (level < max_refinement_levels - 1) {
            refine_contour(&refined_contour, 2);
        }
    }
    
    free(refined_contour.points);
    return j_values;
}

double calculate_path_independence_error(const JIntegralAnalysis *analysis) {
    if (!analysis || analysis->num_results < 2) return 0.0;
    
    // Calculate coefficient of variation (std_dev / mean)
    if (analysis->mean_j_value != 0.0) {
        return analysis->std_deviation / fabs(analysis->mean_j_value);
    }
    
    return 0.0;
}

int validate_j_integral_calculation(const JIntegralAnalysis *analysis, 
                                   double tolerance) {
    if (!analysis) return -1;
    
    double path_error = calculate_path_independence_error(analysis);
    
    if (path_error > tolerance) {
        printf("Warning: Path independence error (%.3f%%) exceeds tolerance (%.3f%%)\n",
               path_error * 100.0, tolerance * 100.0);
        return 1;
    }
    
    // Check for NaN or infinite values
    for (int i = 0; i < analysis->num_results; i++) {
        if (!isfinite(analysis->results[i].j_value)) {
            printf("Error: Invalid J-integral value found\n");
            return -1;
        }
    }
    
    return 0;
}

void export_j_integral_results(const JIntegralAnalysis *analysis,
                              const char *filename, int format) {
    if (!analysis || !filename) return;
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Cannot create output file %s\n", filename);
        return;
    }
    
    if (format == 1) { // CSV format
        fprintf(file, "Contour,J_Value,Energy_Term,Stress_Term,Contour_Length,Num_Points,Computation_Time\n");
        for (int i = 0; i < analysis->num_results; i++) {
            fprintf(file, "%s,%.6e,%.6e,%.6e,%.6f,%d,%.6f\n",
                    analysis->results[i].contour_name,
                    analysis->results[i].j_value,
                    analysis->results[i].energy_term,
                    analysis->results[i].stress_term,
                    analysis->results[i].contour_length,
                    analysis->results[i].num_integration_points,
                    analysis->results[i].computation_time);
        }
    } else { // Text format
        fprintf(file, "J-Integral Analysis Results\n");
        fprintf(file, "===========================\n");
        fprintf(file, "Analysis: %s\n", analysis->analysis_name);
        fprintf(file, "Crack tip: (%.6f, %.6f)\n", analysis->crack_tip.x, analysis->crack_tip.y);
        fprintf(file, "Number of contours: %d\n", analysis->num_results);
        fprintf(file, "Mean J-value: %.6e\n", analysis->mean_j_value);
        fprintf(file, "Standard deviation: %.6e\n", analysis->std_deviation);
        fprintf(file, "Path independence error: %.3f%%\n", analysis->path_independence_error * 100.0);
        fprintf(file, "\nIndividual Results:\n");
        
        for (int i = 0; i < analysis->num_results; i++) {
            fprintf(file, "\nContour: %s\n", analysis->results[i].contour_name);
            fprintf(file, "  J-value: %.6e\n", analysis->results[i].j_value);
            fprintf(file, "  Energy term: %.6e\n", analysis->results[i].energy_term);
            fprintf(file, "  Stress term: %.6e\n", analysis->results[i].stress_term);
            fprintf(file, "  Contour length: %.6f\n", analysis->results[i].contour_length);
            fprintf(file, "  Integration points: %d\n", analysis->results[i].num_integration_points);
            fprintf(file, "  Computation time: %.6f s\n", analysis->results[i].computation_time);
        }
    }
    
    fclose(file);
    printf("J-integral results exported to %s\n", filename);
}

void print_j_integral_result(const JIntegralResult *result) {
    if (!result) return;
    
    printf("\n=== J-INTEGRAL RESULT ===\n");
    printf("Contour: %s\n", result->contour_name);
    printf("J-value: %.6e\n", result->j_value);
    printf("Energy term: %.6e\n", result->energy_term);
    printf("Stress term: %.6e\n", result->stress_term);
    printf("Contour length: %.6f\n", result->contour_length);
    printf("Integration points: %d\n", result->num_integration_points);
    printf("Computation time: %.6f s\n", result->computation_time);
    printf("=========================\n\n");
}

void print_j_integral_analysis(const JIntegralAnalysis *analysis) {
    if (!analysis) return;
    
    printf("\n=== J-INTEGRAL ANALYSIS ===\n");
    printf("Analysis: %s\n", analysis->analysis_name);
    printf("Crack tip: (%.6f, %.6f)\n", analysis->crack_tip.x, analysis->crack_tip.y);
    printf("Number of contours: %d\n", analysis->num_results);
    printf("Mean J-value: %.6e\n", analysis->mean_j_value);
    printf("Standard deviation: %.6e\n", analysis->std_deviation);
    printf("Path independence error: %.3f%%\n", analysis->path_independence_error * 100.0);
    
    printf("\nIndividual Results:\n");
    for (int i = 0; i < analysis->num_results; i++) {
        printf("  %s: J = %.6e\n", analysis->results[i].contour_name, analysis->results[i].j_value);
    }
    printf("===========================\n\n");
}

void calculate_j_integral_statistics(JIntegralAnalysis *analysis) {
    if (!analysis || analysis->num_results == 0) return;
    
    // Calculate mean
    double sum = 0.0;
    for (int i = 0; i < analysis->num_results; i++) {
        sum += analysis->results[i].j_value;
    }
    analysis->mean_j_value = sum / analysis->num_results;
    
    // Calculate standard deviation
    double sum_sq_diff = 0.0;
    for (int i = 0; i < analysis->num_results; i++) {
        double diff = analysis->results[i].j_value - analysis->mean_j_value;
        sum_sq_diff += diff * diff;
    }
    analysis->std_deviation = sqrt(sum_sq_diff / analysis->num_results);
    
    // Calculate path independence error
    analysis->path_independence_error = calculate_path_independence_error(analysis);
}

void free_j_integral_result(JIntegralResult *result) {
    // JIntegralResult contains no dynamically allocated memory
    // This function is provided for consistency
    if (result) {
        memset(result, 0, sizeof(JIntegralResult));
    }
}

void free_j_integral_analysis(JIntegralAnalysis *analysis) {
    if (!analysis) return;
    
    if (analysis->results) {
        free(analysis->results);
    }
    
    free(analysis);
}

int compare_with_analytical_solution(double calculated_j, double analytical_j,
                                    double *relative_error) {
    if (!relative_error) return -1;
    
    if (analytical_j != 0.0) {
        *relative_error = fabs(calculated_j - analytical_j) / fabs(analytical_j);
        return 0;
    }
    
    return -1;
}

double analytical_j_integral_mode_I(double K_I, double E, double nu, int plane_stress) {
    if (K_I <= 0 || E <= 0) return 0.0;
    
    double E_prime;
    if (plane_stress) {
        E_prime = E;
    } else {
        E_prime = E / (1.0 - nu * nu);
    }
    
    return (K_I * K_I) / E_prime;
} 