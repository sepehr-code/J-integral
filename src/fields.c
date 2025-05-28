#define _USE_MATH_DEFINES
#include "../include/fields.h"
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

FieldSet* create_field_set(Mesh *mesh) {
    if (!mesh) {
        fprintf(stderr, "Error: Null mesh pointer\n");
        return NULL;
    }
    
    FieldSet *fields = (FieldSet*)malloc(sizeof(FieldSet));
    if (!fields) {
        fprintf(stderr, "Error: Memory allocation failed for field set\n");
        return NULL;
    }
    
    fields->num_points = mesh->num_nodes;
    fields->mesh = mesh;
    
    fields->data = (FieldData*)malloc(fields->num_points * sizeof(FieldData));
    if (!fields->data) {
        fprintf(stderr, "Error: Memory allocation failed for field data\n");
        free(fields);
        return NULL;
    }
    
    // Initialize field data
    for (int i = 0; i < fields->num_points; i++) {
        fields->data[i].node_id = i;
        fields->data[i].coord = mesh->nodes[i].coord;
        
        // Initialize all fields to zero
        fields->data[i].disp.ux = 0.0;
        fields->data[i].disp.uy = 0.0;
        fields->data[i].stress.xx = 0.0;
        fields->data[i].stress.yy = 0.0;
        fields->data[i].stress.xy = 0.0;
        fields->data[i].strain.xx = 0.0;
        fields->data[i].strain.yy = 0.0;
        fields->data[i].strain.xy = 0.0;
        fields->data[i].strain_energy_density = 0.0;
    }
    
    return fields;
}

void generate_synthetic_stress_field(FieldSet *fields, Point2D crack_tip, 
                                   double applied_stress) {
    if (!fields) return;
    
    printf("Generating synthetic stress field...\n");
    
    for (int i = 0; i < fields->num_points; i++) {
        double x = fields->data[i].coord.x - crack_tip.x;
        double y = fields->data[i].coord.y - crack_tip.y;
        double r = sqrt(x * x + y * y);
        
        if (r < 1e-10) r = 1e-10; // Avoid singularity at crack tip
        
        // Simple polynomial stress field for testing
        // This is not physically accurate but useful for validation
        fields->data[i].stress.xx = applied_stress * (1.0 + 0.1 * x + 0.05 * y);
        fields->data[i].stress.yy = applied_stress * (0.8 + 0.05 * x + 0.1 * y);
        fields->data[i].stress.xy = applied_stress * (0.1 * x * y / (r + 1.0));
        
        // Simple displacement field (assuming linear elastic)
        double E = 200000.0; // Young's modulus (MPa)
        double nu = 0.3;     // Poisson's ratio
        
        fields->data[i].disp.ux = (fields->data[i].stress.xx - nu * fields->data[i].stress.yy) * x / E;
        fields->data[i].disp.uy = (fields->data[i].stress.yy - nu * fields->data[i].stress.xx) * y / E;
        
        // Calculate strain from stress (plane stress assumption)
        fields->data[i].strain.xx = (fields->data[i].stress.xx - nu * fields->data[i].stress.yy) / E;
        fields->data[i].strain.yy = (fields->data[i].stress.yy - nu * fields->data[i].stress.xx) / E;
        fields->data[i].strain.xy = fields->data[i].stress.xy / (E / (2.0 * (1.0 + nu)));
        
        // Calculate strain energy density
        fields->data[i].strain_energy_density = 0.5 * (
            fields->data[i].stress.xx * fields->data[i].strain.xx +
            fields->data[i].stress.yy * fields->data[i].strain.yy +
            2.0 * fields->data[i].stress.xy * fields->data[i].strain.xy
        );
    }
    
    printf("Synthetic stress field generated for %d points\n", fields->num_points);
}

void generate_williams_field(FieldSet *fields, Point2D crack_tip, 
                           double K_I, int num_terms) {
    if (!fields) return;
    
    printf("Generating Williams series field (K_I = %.2f, %d terms)...\n", K_I, num_terms);
    
    // Material properties (typical steel)
    double E = 200000.0; // Young's modulus (MPa)
    double nu = 0.3;     // Poisson's ratio
    double mu = E / (2.0 * (1.0 + nu)); // Shear modulus
    double kappa = 3.0 - 4.0 * nu; // Plane strain condition
    
    for (int i = 0; i < fields->num_points; i++) {
        double x = fields->data[i].coord.x - crack_tip.x;
        double y = fields->data[i].coord.y - crack_tip.y;
        double r = sqrt(x * x + y * y);
        double theta = atan2(y, x);
        
        if (r < 1e-10) r = 1e-10; // Avoid singularity
        
        // Mode I Williams series (first term only for simplicity)
        double sqrt_r = sqrt(r);
        double cos_half_theta = cos(theta / 2.0);
        double sin_half_theta = sin(theta / 2.0);
        double cos_3half_theta = cos(3.0 * theta / 2.0);
        double sin_3half_theta = sin(3.0 * theta / 2.0);
        
        // Stress field (Mode I, first term)
        double coeff = K_I / sqrt(2.0 * M_PI * r);
        
        fields->data[i].stress.xx = coeff * cos_half_theta * (1.0 - sin_half_theta * sin_3half_theta);
        fields->data[i].stress.yy = coeff * cos_half_theta * (1.0 + sin_half_theta * sin_3half_theta);
        fields->data[i].stress.xy = coeff * sin_half_theta * cos_half_theta * cos_3half_theta;
        
        // Displacement field (Mode I, first term)
        double disp_coeff = K_I / (2.0 * mu) * sqrt(r / (2.0 * M_PI));
        
        fields->data[i].disp.ux = disp_coeff * cos_half_theta * (kappa - 1.0 + 2.0 * sin_half_theta * sin_half_theta);
        fields->data[i].disp.uy = disp_coeff * sin_half_theta * (kappa + 1.0 - 2.0 * cos_half_theta * cos_half_theta);
        
        // Calculate strain from stress
        fields->data[i].strain.xx = (fields->data[i].stress.xx - nu * fields->data[i].stress.yy) / E;
        fields->data[i].strain.yy = (fields->data[i].stress.yy - nu * fields->data[i].stress.xx) / E;
        fields->data[i].strain.xy = fields->data[i].stress.xy / (2.0 * mu);
        
        // Calculate strain energy density
        fields->data[i].strain_energy_density = 0.5 * (
            fields->data[i].stress.xx * fields->data[i].strain.xx +
            fields->data[i].stress.yy * fields->data[i].strain.yy +
            2.0 * fields->data[i].stress.xy * fields->data[i].strain.xy
        );
    }
    
    printf("Williams series field generated for %d points\n", fields->num_points);
}

FieldSet* load_fields_from_file(const char *filename, Mesh *mesh) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open field file %s\n", filename);
        return NULL;
    }
    
    FieldSet *fields = create_field_set(mesh);
    if (!fields) {
        fclose(file);
        return NULL;
    }
    
    char line[256];
    int point_count = 0;
    
    // Skip header lines starting with #
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#') {
            // Parse field data line
            int node_id;
            double x, y, ux, uy, sxx, syy, sxy;
            
            if (sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf", 
                      &node_id, &x, &y, &ux, &uy, &sxx, &syy, &sxy) == 8) {
                
                if (point_count < fields->num_points) {
                    fields->data[point_count].node_id = node_id;
                    fields->data[point_count].coord.x = x;
                    fields->data[point_count].coord.y = y;
                    fields->data[point_count].disp.ux = ux;
                    fields->data[point_count].disp.uy = uy;
                    fields->data[point_count].stress.xx = sxx;
                    fields->data[point_count].stress.yy = syy;
                    fields->data[point_count].stress.xy = sxy;
                    
                    point_count++;
                }
            }
        }
    }
    
    fclose(file);
    
    // Update strain energy density
    update_strain_energy_density(fields);
    
    printf("Loaded field data for %d points from %s\n", point_count, filename);
    return fields;
}

void save_fields_to_file(const FieldSet *fields, const char *filename) {
    if (!fields || !filename) return;
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Cannot create field file %s\n", filename);
        return;
    }
    
    // Write header
    fprintf(file, "# Field data file\n");
    fprintf(file, "# node_id x y ux uy sigma_xx sigma_yy sigma_xy strain_energy_density\n");
    
    // Write field data
    for (int i = 0; i < fields->num_points; i++) {
        fprintf(file, "%d %.6f %.6f %.6e %.6e %.6e %.6e %.6e %.6e\n",
                fields->data[i].node_id,
                fields->data[i].coord.x, fields->data[i].coord.y,
                fields->data[i].disp.ux, fields->data[i].disp.uy,
                fields->data[i].stress.xx, fields->data[i].stress.yy, fields->data[i].stress.xy,
                fields->data[i].strain_energy_density);
    }
    
    fclose(file);
    printf("Field data saved to %s\n", filename);
}

int interpolate_fields(const FieldSet *fields, Point2D point,
                      Displacement *disp, Stress *stress, 
                      double *strain_energy_density) {
    if (!fields || !disp || !stress || !strain_energy_density) return -1;
    
    // Find the element containing the point
    int elem_id = find_element_containing_point(fields->mesh, point);
    if (elem_id < 0) return -1;
    
    // Get element nodes
    Element *elem = &fields->mesh->elements[elem_id];
    
    // Simple bilinear interpolation using element corner nodes
    double dx = fields->mesh->width / (fields->mesh->nx - 1);
    double dy = fields->mesh->height / (fields->mesh->ny - 1);
    
    // Local coordinates within element
    int i = elem_id % (fields->mesh->nx - 1);
    int j = elem_id / (fields->mesh->nx - 1);
    
    double x0 = i * dx;
    double y0 = j * dy;
    double xi = (point.x - x0) / dx;
    double eta = (point.y - y0) / dy;
    
    // Bilinear shape functions
    double N[4];
    N[0] = (1.0 - xi) * (1.0 - eta);  // bottom-left
    N[1] = xi * (1.0 - eta);          // bottom-right
    N[2] = xi * eta;                  // top-right
    N[3] = (1.0 - xi) * eta;          // top-left
    
    // Initialize interpolated values
    disp->ux = 0.0; disp->uy = 0.0;
    stress->xx = 0.0; stress->yy = 0.0; stress->xy = 0.0;
    *strain_energy_density = 0.0;
    
    // Interpolate field values
    for (int k = 0; k < 4; k++) {
        int node_id = elem->nodes[k];
        if (node_id < fields->num_points) {
            disp->ux += N[k] * fields->data[node_id].disp.ux;
            disp->uy += N[k] * fields->data[node_id].disp.uy;
            stress->xx += N[k] * fields->data[node_id].stress.xx;
            stress->yy += N[k] * fields->data[node_id].stress.yy;
            stress->xy += N[k] * fields->data[node_id].stress.xy;
            *strain_energy_density += N[k] * fields->data[node_id].strain_energy_density;
        }
    }
    
    return 0;
}

Strain calculate_strain_from_gradients(double du_dx, double du_dy,
                                     double dv_dx, double dv_dy) {
    Strain strain;
    strain.xx = du_dx;
    strain.yy = dv_dy;
    strain.xy = 0.5 * (du_dy + dv_dx);
    return strain;
}

double calculate_strain_energy_density(const Stress *stress, const Strain *strain) {
    if (!stress || !strain) return 0.0;
    
    return 0.5 * (stress->xx * strain->xx + 
                  stress->yy * strain->yy + 
                  2.0 * stress->xy * strain->xy);
}

int calculate_displacement_gradients(const FieldSet *fields, Point2D point,
                                   double *du_dx, double *du_dy,
                                   double *dv_dx, double *dv_dy) {
    if (!fields || !du_dx || !du_dy || !dv_dx || !dv_dy) return -1;
    
    // Simple finite difference approximation
    double h = 1e-6; // Small perturbation
    
    Point2D p_plus_x = {point.x + h, point.y};
    Point2D p_minus_x = {point.x - h, point.y};
    Point2D p_plus_y = {point.x, point.y + h};
    Point2D p_minus_y = {point.x, point.y - h};
    
    Displacement disp_plus_x, disp_minus_x, disp_plus_y, disp_minus_y;
    Stress stress_dummy;
    double sed_dummy;
    
    // Calculate gradients using central differences
    if (interpolate_fields(fields, p_plus_x, &disp_plus_x, &stress_dummy, &sed_dummy) != 0 ||
        interpolate_fields(fields, p_minus_x, &disp_minus_x, &stress_dummy, &sed_dummy) != 0 ||
        interpolate_fields(fields, p_plus_y, &disp_plus_y, &stress_dummy, &sed_dummy) != 0 ||
        interpolate_fields(fields, p_minus_y, &disp_minus_y, &stress_dummy, &sed_dummy) != 0) {
        return -1;
    }
    
    *du_dx = (disp_plus_x.ux - disp_minus_x.ux) / (2.0 * h);
    *du_dy = (disp_plus_y.ux - disp_minus_y.ux) / (2.0 * h);
    *dv_dx = (disp_plus_x.uy - disp_minus_x.uy) / (2.0 * h);
    *dv_dy = (disp_plus_y.uy - disp_minus_y.uy) / (2.0 * h);
    
    return 0;
}

void update_strain_energy_density(FieldSet *fields) {
    if (!fields) return;
    
    for (int i = 0; i < fields->num_points; i++) {
        fields->data[i].strain_energy_density = calculate_strain_energy_density(
            &fields->data[i].stress, &fields->data[i].strain);
    }
}

void free_field_set(FieldSet *fields) {
    if (!fields) return;
    
    if (fields->data) {
        free(fields->data);
    }
    
    free(fields);
}

void print_field_info(const FieldSet *fields, int point_id) {
    if (!fields || point_id < 0 || point_id >= fields->num_points) {
        printf("Error: Invalid field set or point ID\n");
        return;
    }
    
    FieldData *data = &fields->data[point_id];
    
    printf("\n=== FIELD DATA AT POINT %d ===\n", point_id);
    printf("Coordinates: (%.6f, %.6f)\n", data->coord.x, data->coord.y);
    printf("Displacement: ux=%.6e, uy=%.6e\n", data->disp.ux, data->disp.uy);
    printf("Stress: σxx=%.6e, σyy=%.6e, σxy=%.6e\n", 
           data->stress.xx, data->stress.yy, data->stress.xy);
    printf("Strain: εxx=%.6e, εyy=%.6e, εxy=%.6e\n", 
           data->strain.xx, data->strain.yy, data->strain.xy);
    printf("Strain Energy Density: %.6e\n", data->strain_energy_density);
    printf("==============================\n\n");
}

int validate_field_data(const FieldSet *fields) {
    if (!fields || !fields->data) return -1;
    
    int errors = 0;
    
    for (int i = 0; i < fields->num_points; i++) {
        // Check for NaN or infinite values
        if (!isfinite(fields->data[i].disp.ux) || !isfinite(fields->data[i].disp.uy) ||
            !isfinite(fields->data[i].stress.xx) || !isfinite(fields->data[i].stress.yy) ||
            !isfinite(fields->data[i].stress.xy) || !isfinite(fields->data[i].strain_energy_density)) {
            errors++;
        }
    }
    
    if (errors > 0) {
        printf("Warning: Found %d points with invalid field values\n", errors);
    }
    
    return errors;
} 