#include "../include/j_integral.h"
#include <string.h>

void print_usage(const char *program_name) {
    printf("Usage: %s [options] [mesh_file] [field_file]\n", program_name);
    printf("\nOptions:\n");
    printf("  -h, --help          Show this help message\n");
    printf("  -s, --synthetic     Use synthetic field data (default)\n");
    printf("  -w, --williams      Use Williams series field\n");
    printf("  -K <value>          Stress intensity factor for Williams series (default: 1000)\n");
    printf("  -r <radius>         Contour radius (default: 0.1)\n");
    printf("  -n <points>         Number of contour points (default: 64)\n");
    printf("  -c <num>            Number of contours for path independence study (default: 3)\n");
    printf("  -o <filename>       Output file for results\n");
    printf("  --csv               Export results in CSV format\n");
    printf("  --convergence       Perform convergence study\n");
    printf("  --validate          Compare with analytical solution\n");
    printf("\nExamples:\n");
    printf("  %s                           # Run with default synthetic data\n", program_name);
    printf("  %s -w -K 1500 -r 0.05        # Williams series with K_I=1500, radius=0.05\n", program_name);
    printf("  %s mesh.dat field.dat        # Load data from files\n", program_name);
    printf("  %s --convergence --validate  # Full analysis with validation\n", program_name);
}

void run_basic_example() {
    printf("\n=== BASIC J-INTEGRAL EXAMPLE ===\n");
    
    // Create a simple mesh
    Mesh *mesh = create_rectangular_mesh(21, 21, 2.0, 2.0);
    if (!mesh) {
        fprintf(stderr, "Error: Failed to create mesh\n");
        return;
    }
    
    // Define a crack at the center
    Point2D crack_tip = {1.0, 1.0};
    define_crack(mesh, 0.5, 1.0, 1.0, 1.0);
    
    // Create field set
    FieldSet *fields = create_field_set(mesh);
    if (!fields) {
        fprintf(stderr, "Error: Failed to create field set\n");
        free_mesh(mesh);
        return;
    }
    
    // Generate synthetic stress field
    generate_synthetic_stress_field(fields, crack_tip, 100.0);
    
    // Create contour set with multiple contours for path independence study
    ContourSet *contour_set = create_contour_set(crack_tip, 3);
    add_circular_contour(contour_set, 0.05, 32, "Inner_Circle");
    add_circular_contour(contour_set, 0.1, 64, "Middle_Circle");
    add_circular_contour(contour_set, 0.15, 96, "Outer_Circle");
    
    // Calculate J-integral for multiple contours
    JIntegralAnalysis *analysis = calculate_j_integral_multiple_contours(fields, contour_set,
                                                                        INTEGRATION_MIDPOINT,
                                                                        "Basic_Example");
    
    // Print results
    print_j_integral_analysis(analysis);
    
    // Validate path independence
    double tolerance = 0.05; // 5% tolerance
    int validation_result = validate_j_integral_calculation(analysis, tolerance);
    if (validation_result == 0) {
        printf("✓ Path independence validated (error < %.1f%%)\n", tolerance * 100.0);
    } else if (validation_result == 1) {
        printf("⚠ Path independence warning (error > %.1f%%)\n", tolerance * 100.0);
    } else {
        printf("✗ Validation failed\n");
    }
    
    // Cleanup
    free_j_integral_analysis(analysis);
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    
    printf("=== BASIC EXAMPLE COMPLETED ===\n\n");
}

void run_williams_example(double K_I, double contour_radius, int num_points) {
    printf("\n=== WILLIAMS SERIES EXAMPLE ===\n");
    printf("K_I = %.1f MPa√m, radius = %.3f, points = %d\n", K_I, contour_radius, num_points);
    
    // Create mesh
    Mesh *mesh = create_rectangular_mesh(41, 41, 1.0, 1.0);
    if (!mesh) {
        fprintf(stderr, "Error: Failed to create mesh\n");
        return;
    }
    
    // Define crack
    Point2D crack_tip = {0.5, 0.5};
    define_crack(mesh, 0.2, 0.5, 0.5, 0.5);
    
    // Create field set
    FieldSet *fields = create_field_set(mesh);
    if (!fields) {
        fprintf(stderr, "Error: Failed to create field set\n");
        free_mesh(mesh);
        return;
    }
    
    // Generate Williams series field
    generate_williams_field(fields, crack_tip, K_I, 1);
    
    // Create single contour
    ContourSet *contour_set = create_contour_set(crack_tip, 1);
    add_circular_contour(contour_set, contour_radius, num_points, "Williams_Contour");
    
    // Calculate J-integral
    JIntegralResult result = calculate_j_integral_single_contour(fields,
                                                               &contour_set->contours[0],
                                                               INTEGRATION_MIDPOINT);
    
    // Print result
    print_j_integral_result(&result);
    
    // Compare with analytical solution
    double E = 200000.0; // MPa
    double nu = 0.3;
    double analytical_j = analytical_j_integral_mode_I(K_I, E, nu, 0); // plane strain
    double relative_error;
    
    if (compare_with_analytical_solution(result.j_value, analytical_j, &relative_error) == 0) {
        printf("Analytical J-integral: %.6e\n", analytical_j);
        printf("Relative error: %.3f%%\n", relative_error * 100.0);
        
        if (relative_error < 0.1) {
            printf("✓ Excellent agreement with analytical solution\n");
        } else if (relative_error < 0.2) {
            printf("✓ Good agreement with analytical solution\n");
        } else {
            printf("⚠ Moderate agreement with analytical solution\n");
        }
    }
    
    // Cleanup
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    
    printf("=== WILLIAMS EXAMPLE COMPLETED ===\n\n");
}

void run_convergence_study() {
    printf("\n=== CONVERGENCE STUDY ===\n");
    
    // Create base mesh
    Mesh *mesh = create_rectangular_mesh(11, 11, 1.0, 1.0);
    if (!mesh) {
        fprintf(stderr, "Error: Failed to create base mesh\n");
        return;
    }
    
    Point2D crack_tip = {0.5, 0.5};
    define_crack(mesh, 0.2, 0.5, 0.5, 0.5);
    
    // Create field set
    FieldSet *fields = create_field_set(mesh);
    if (!fields) {
        fprintf(stderr, "Error: Failed to create field set\n");
        free_mesh(mesh);
        return;
    }
    
    // Generate Williams series field for accurate convergence study
    generate_williams_field(fields, crack_tip, 1000.0, 1);
    
    // Mesh refinement convergence study
    printf("\n--- Mesh Refinement Convergence ---\n");
    double *mesh_j_values = convergence_study_mesh_refinement(mesh, fields, crack_tip, 4, 0.1);
    if (mesh_j_values) {
        printf("Mesh convergence completed\n");
        free(mesh_j_values);
    }
    
    // Contour refinement convergence study
    printf("\n--- Contour Refinement Convergence ---\n");
    ContourSet *contour_set = create_contour_set(crack_tip, 1);
    add_circular_contour(contour_set, 0.1, 16, "Base_Contour");
    
    double *contour_j_values = convergence_study_contour_refinement(fields,
                                                                   &contour_set->contours[0], 4);
    if (contour_j_values) {
        printf("Contour convergence completed\n");
        free(contour_j_values);
    }
    
    // Cleanup
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    
    printf("=== CONVERGENCE STUDY COMPLETED ===\n\n");
}

void run_three_point_bend_example() {
    printf("\n=== THREE-POINT BEND TEST SIMULATION ===\n");
    
    // Create specimen geometry (typical 3-point bend specimen)
    double specimen_length = 0.1;  // 100 mm
    double specimen_height = 0.08; // 80 mm
    double crack_length = 0.01;    // 10 mm (a/W = 0.5)
    
    Mesh *mesh = create_rectangular_mesh(51, 21, specimen_length, specimen_height);
    if (!mesh) {
        fprintf(stderr, "Error: Failed to create specimen mesh\n");
        return;
    }
    
    // Define edge crack (typical for 3-point bend)
    Point2D crack_tip = {crack_length, specimen_height / 2.0};
    define_crack(mesh, 0.0, specimen_height / 2.0, crack_length, specimen_height / 2.0);
    
    // Create field set
    FieldSet *fields = create_field_set(mesh);
    if (!fields) {
        fprintf(stderr, "Error: Failed to create field set\n");
        free_mesh(mesh);
        return;
    }
    
    // Generate stress field (simplified 3-point bend loading)
    // In reality, this would come from FEM analysis
    double applied_load = 1000.0; // N
    double moment = applied_load * specimen_length / 4.0; // Maximum moment at center
    double I = specimen_height * specimen_height * specimen_height / 12.0; // Second moment of area
    
    for (int i = 0; i < fields->num_points; i++) {
        double x = fields->data[i].coord.x;
        double y = fields->data[i].coord.y - specimen_height / 2.0; // Center at y=0
        
        // Beam bending stress (simplified)
        double sigma_xx = moment * y / I;
        double sigma_yy = 0.0;
        double sigma_xy = 0.0;
        
        // Add stress concentration near crack tip
        double r = point_distance(fields->data[i].coord, crack_tip);
        if (r > 1e-6) {
            double stress_concentration = 1.0 + 2.0 * exp(-r / 0.002); // Exponential decay
            sigma_xx *= stress_concentration;
        }
        
        fields->data[i].stress.xx = sigma_xx;
        fields->data[i].stress.yy = sigma_yy;
        fields->data[i].stress.xy = sigma_xy;
        
        // Calculate strain (plane stress)
        double E = 200000.0; // MPa
        double nu = 0.3;
        fields->data[i].strain.xx = (sigma_xx - nu * sigma_yy) / E;
        fields->data[i].strain.yy = (sigma_yy - nu * sigma_xx) / E;
        fields->data[i].strain.xy = sigma_xy / (E / (2.0 * (1.0 + nu)));
        
        // Strain energy density
        fields->data[i].strain_energy_density = 0.5 * (
            sigma_xx * fields->data[i].strain.xx +
            sigma_yy * fields->data[i].strain.yy +
            2.0 * sigma_xy * fields->data[i].strain.xy
        );
    }
    
    // Create multiple contours around crack tip
    ContourSet *contour_set = create_contour_set(crack_tip, 4);
    add_circular_contour(contour_set, 0.001, 32, "r=1mm");
    add_circular_contour(contour_set, 0.002, 48, "r=2mm");
    add_circular_contour(contour_set, 0.003, 64, "r=3mm");
    add_rectangular_contour(contour_set, 0.008, 0.006, 20, "Rectangle");
    
    // Calculate J-integral
    JIntegralAnalysis *analysis = calculate_j_integral_multiple_contours(fields, contour_set,
                                                                        INTEGRATION_MIDPOINT,
                                                                        "Three_Point_Bend");
    
    // Print results
    print_j_integral_analysis(analysis);
    
    // Calculate fracture toughness (K_IC) from J-integral
    double E = 200000.0; // MPa
    double nu = 0.3;
    double E_prime = E / (1.0 - nu * nu); // Plane strain
    double K_IC = sqrt(analysis->mean_j_value * E_prime);
    
    printf("Calculated fracture toughness K_IC: %.1f MPa√m\n", K_IC);
    
    // Export results
    export_j_integral_results(analysis, "three_point_bend_results.txt", 0);
    export_j_integral_results(analysis, "three_point_bend_results.csv", 1);
    
    // Cleanup
    free_j_integral_analysis(analysis);
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    
    printf("=== THREE-POINT BEND SIMULATION COMPLETED ===\n\n");
}

int main(int argc, char *argv[]) {
    printf("J-Integral Calculation for Fracture Mechanics\n");
    printf("==============================================\n");
    printf("PhD Research Implementation (2010)\n\n");
    
    // Default parameters
    int use_synthetic = 1;
    int use_williams = 0;
    double K_I = 1000.0;
    double contour_radius = 0.1;
    int num_points = 64;
    int num_contours = 3;
    char *output_file = NULL;
    int csv_format = 0;
    int run_convergence = 0;
    int run_validation = 0;
    char *mesh_file = NULL;
    char *field_file = NULL;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--synthetic") == 0) {
            use_synthetic = 1;
            use_williams = 0;
        } else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--williams") == 0) {
            use_williams = 1;
            use_synthetic = 0;
        } else if (strcmp(argv[i], "-K") == 0 && i + 1 < argc) {
            K_I = atof(argv[++i]);
        } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
            contour_radius = atof(argv[++i]);
        } else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            num_points = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            num_contours = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--csv") == 0) {
            csv_format = 1;
        } else if (strcmp(argv[i], "--convergence") == 0) {
            run_convergence = 1;
        } else if (strcmp(argv[i], "--validate") == 0) {
            run_validation = 1;
        } else if (argv[i][0] != '-') {
            // Non-option arguments are treated as file names
            if (!mesh_file) {
                mesh_file = argv[i];
            } else if (!field_file) {
                field_file = argv[i];
            }
        }
    }
    
    // Run examples based on command line options
    if (run_convergence) {
        run_convergence_study();
    }
    
    if (use_williams) {
        run_williams_example(K_I, contour_radius, num_points);
        
        if (run_validation) {
            printf("Williams series validation completed above.\n\n");
        }
    } else if (mesh_file && field_file) {
        printf("Loading data from files: %s, %s\n", mesh_file, field_file);
        
        // Load mesh and fields from files
        Mesh *mesh = load_mesh_from_file(mesh_file);
        if (!mesh) {
            fprintf(stderr, "Error: Failed to load mesh from %s\n", mesh_file);
            return 1;
        }
        
        FieldSet *fields = load_fields_from_file(field_file, mesh);
        if (!fields) {
            fprintf(stderr, "Error: Failed to load fields from %s\n", field_file);
            free_mesh(mesh);
            return 1;
        }
        
        // Use crack tip from mesh
        Point2D crack_tip = mesh->crack.tip;
        
        // Create contour set
        ContourSet *contour_set = create_contour_set(crack_tip, num_contours);
        for (int i = 0; i < num_contours; i++) {
            char contour_name[64];
            snprintf(contour_name, sizeof(contour_name), "Contour_%d", i + 1);
            double radius = contour_radius * (1.0 + 0.5 * i);
            add_circular_contour(contour_set, radius, num_points, contour_name);
        }
        
        // Calculate J-integral
        JIntegralAnalysis *analysis = calculate_j_integral_multiple_contours(fields, contour_set,
                                                                            INTEGRATION_MIDPOINT,
                                                                            "File_Data_Analysis");
        
        // Print results
        print_j_integral_analysis(analysis);
        
        // Export if requested
        if (output_file) {
            export_j_integral_results(analysis, output_file, csv_format);
        }
        
        // Cleanup
        free_j_integral_analysis(analysis);
        free_contour_set(contour_set);
        free_field_set(fields);
        free_mesh(mesh);
        
    } else {
        // Run default examples
        run_basic_example();
        run_three_point_bend_example();
    }
    
    printf("J-Integral calculation completed successfully.\n");
    printf("For more options, run: %s --help\n", argv[0]);
    
    return 0;
} 