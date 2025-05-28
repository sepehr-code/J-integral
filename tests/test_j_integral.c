#include "../include/j_integral.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define TEST_TOLERANCE 1e-6
#define ASSERT_NEAR(a, b, tol) assert(fabs((a) - (b)) < (tol))

void test_mesh_creation() {
    printf("Testing mesh creation...\n");
    
    Mesh *mesh = create_rectangular_mesh(11, 11, 1.0, 1.0);
    assert(mesh != NULL);
    assert(mesh->nx == 11);
    assert(mesh->ny == 11);
    assert(mesh->num_nodes == 121);
    assert(mesh->num_elements == 100);
    assert(mesh->width == 1.0);
    assert(mesh->height == 1.0);
    
    free_mesh(mesh);
    printf("✓ Mesh creation test passed\n");
}

void test_field_creation() {
    printf("Testing field creation...\n");
    
    Mesh *mesh = create_rectangular_mesh(5, 5, 1.0, 1.0);
    assert(mesh != NULL);
    
    FieldSet *fields = create_field_set(mesh);
    assert(fields != NULL);
    assert(fields->num_points == 25);
    assert(fields->mesh == mesh);
    
    free_field_set(fields);
    free_mesh(mesh);
    printf("✓ Field creation test passed\n");
}

void test_contour_creation() {
    printf("Testing contour creation...\n");
    
    Point2D crack_tip = {0.5, 0.5};
    ContourSet *contour_set = create_contour_set(crack_tip, 2);
    assert(contour_set != NULL);
    
    int result = add_circular_contour(contour_set, 0.1, 32, "test_circle");
    assert(result == 0);
    assert(contour_set->num_contours == 1);
    
    result = add_rectangular_contour(contour_set, 0.2, 0.2, 10, "test_rectangle");
    assert(result == 0);
    assert(contour_set->num_contours == 2);
    
    free_contour_set(contour_set);
    printf("✓ Contour creation test passed\n");
}

void test_williams_field_validation() {
    printf("Testing Williams series field validation...\n");
    
    // Create small mesh
    Mesh *mesh = create_rectangular_mesh(21, 21, 1.0, 1.0);
    Point2D crack_tip = {0.5, 0.5};
    define_crack(mesh, 0.2, 0.5, 0.5, 0.5);
    
    FieldSet *fields = create_field_set(mesh);
    
    // Generate Williams field with known K_I
    double K_I = 1000.0; // MPa√m
    generate_williams_field(fields, crack_tip, K_I, 1);
    
    // Create contour
    ContourSet *contour_set = create_contour_set(crack_tip, 1);
    add_circular_contour(contour_set, 0.1, 64, "validation_contour");
    
    // Calculate J-integral
    JIntegralResult result = calculate_j_integral_single_contour(fields,
                                                               &contour_set->contours[0],
                                                               INTEGRATION_MIDPOINT);
    
    // Compare with analytical solution
    double E = 200000.0; // MPa
    double nu = 0.3;
    double analytical_j = analytical_j_integral_mode_I(K_I, E, nu, 0); // plane strain
    
    double relative_error = fabs(result.j_value - analytical_j) / analytical_j;
    printf("  Calculated J: %.6e\n", result.j_value);
    printf("  Analytical J: %.6e\n", analytical_j);
    printf("  Relative error: %.3f%%\n", relative_error * 100.0);
    
    // Should be within 20% for this simple test
    assert(relative_error < 0.2);
    
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    printf("✓ Williams field validation test passed\n");
}

void test_path_independence() {
    printf("Testing path independence...\n");
    
    // Create mesh and fields
    Mesh *mesh = create_rectangular_mesh(31, 31, 2.0, 2.0);
    Point2D crack_tip = {1.0, 1.0};
    define_crack(mesh, 0.5, 1.0, 1.0, 1.0);
    
    FieldSet *fields = create_field_set(mesh);
    generate_williams_field(fields, crack_tip, 1000.0, 1);
    
    // Create multiple contours
    ContourSet *contour_set = create_contour_set(crack_tip, 3);
    add_circular_contour(contour_set, 0.05, 32, "inner");
    add_circular_contour(contour_set, 0.1, 64, "middle");
    add_circular_contour(contour_set, 0.15, 96, "outer");
    
    // Calculate J-integral for all contours
    JIntegralAnalysis *analysis = calculate_j_integral_multiple_contours(fields, contour_set,
                                                                        INTEGRATION_MIDPOINT,
                                                                        "path_independence_test");
    
    // Check path independence (should be within 10% for Williams field)
    double path_error = calculate_path_independence_error(analysis);
    printf("  Path independence error: %.3f%%\n", path_error * 100.0);
    
    assert(path_error < 0.1); // 10% tolerance
    
    free_j_integral_analysis(analysis);
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    printf("✓ Path independence test passed\n");
}

void test_analytical_j_integral() {
    printf("Testing analytical J-integral formula...\n");
    
    double K_I = 1000.0; // MPa√m
    double E = 200000.0; // MPa
    double nu = 0.3;
    
    // Plane stress
    double J_plane_stress = analytical_j_integral_mode_I(K_I, E, nu, 1);
    double expected_plane_stress = (K_I * K_I) / E;
    ASSERT_NEAR(J_plane_stress, expected_plane_stress, TEST_TOLERANCE);
    
    // Plane strain
    double J_plane_strain = analytical_j_integral_mode_I(K_I, E, nu, 0);
    double E_prime = E / (1.0 - nu * nu);
    double expected_plane_strain = (K_I * K_I) / E_prime;
    ASSERT_NEAR(J_plane_strain, expected_plane_strain, TEST_TOLERANCE);
    
    printf("✓ Analytical J-integral formula test passed\n");
}

void test_integration_methods() {
    printf("Testing different integration methods...\n");
    
    // Create simple test case
    Mesh *mesh = create_rectangular_mesh(21, 21, 1.0, 1.0);
    Point2D crack_tip = {0.5, 0.5};
    FieldSet *fields = create_field_set(mesh);
    generate_synthetic_stress_field(fields, crack_tip, 100.0);
    
    ContourSet *contour_set = create_contour_set(crack_tip, 1);
    add_circular_contour(contour_set, 0.1, 64, "test_contour");
    
    // Test all integration methods
    double j_midpoint = j_integral_midpoint_rule(fields, &contour_set->contours[0]);
    double j_trapezoidal = j_integral_trapezoidal_rule(fields, &contour_set->contours[0]);
    double j_simpson = j_integral_simpson_rule(fields, &contour_set->contours[0]);
    
    printf("  Midpoint: %.6e\n", j_midpoint);
    printf("  Trapezoidal: %.6e\n", j_trapezoidal);
    printf("  Simpson: %.6e\n", j_simpson);
    
    // All methods should give finite results
    assert(isfinite(j_midpoint));
    assert(isfinite(j_trapezoidal));
    assert(isfinite(j_simpson));
    
    // Results should be reasonably close (within order of magnitude)
    assert(fabs(j_trapezoidal - j_midpoint) / fabs(j_midpoint) < 1.0);
    assert(fabs(j_simpson - j_midpoint) / fabs(j_midpoint) < 1.0);
    
    free_contour_set(contour_set);
    free_field_set(fields);
    free_mesh(mesh);
    printf("✓ Integration methods test passed\n");
}

void test_memory_management() {
    printf("Testing memory management...\n");
    
    // Test multiple allocations and deallocations
    for (int i = 0; i < 10; i++) {
        Mesh *mesh = create_rectangular_mesh(10, 10, 1.0, 1.0);
        FieldSet *fields = create_field_set(mesh);
        ContourSet *contour_set = create_contour_set((Point2D){0.5, 0.5}, 2);
        
        add_circular_contour(contour_set, 0.1, 32, "test");
        add_rectangular_contour(contour_set, 0.2, 0.2, 10, "test2");
        
        free_contour_set(contour_set);
        free_field_set(fields);
        free_mesh(mesh);
    }
    
    printf("✓ Memory management test passed\n");
}

int main() {
    printf("J-Integral Test Suite\n");
    printf("====================\n\n");
    
    test_mesh_creation();
    test_field_creation();
    test_contour_creation();
    test_analytical_j_integral();
    test_integration_methods();
    test_williams_field_validation();
    test_path_independence();
    test_memory_management();
    
    printf("\n=== ALL TESTS PASSED ===\n");
    printf("J-integral implementation validated successfully!\n");
    
    return 0;
} 