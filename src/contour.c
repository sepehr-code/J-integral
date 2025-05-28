#define _USE_MATH_DEFINES
#include "../include/contour.h"
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ContourSet* create_contour_set(Point2D crack_tip, int num_contours) {
    if (num_contours <= 0) {
        fprintf(stderr, "Error: Invalid number of contours\n");
        return NULL;
    }
    
    ContourSet *contour_set = (ContourSet*)malloc(sizeof(ContourSet));
    if (!contour_set) {
        fprintf(stderr, "Error: Memory allocation failed for contour set\n");
        return NULL;
    }
    
    contour_set->num_contours = 0; // Will be incremented as contours are added
    contour_set->crack_tip = crack_tip;
    
    contour_set->contours = (Contour*)malloc(num_contours * sizeof(Contour));
    if (!contour_set->contours) {
        fprintf(stderr, "Error: Memory allocation failed for contours\n");
        free(contour_set);
        return NULL;
    }
    
    return contour_set;
}

int add_circular_contour(ContourSet *contour_set, double radius, 
                        int num_points, const char *name) {
    if (!contour_set || radius <= 0 || num_points < 3) {
        fprintf(stderr, "Error: Invalid parameters for circular contour\n");
        return -1;
    }
    
    int contour_id = contour_set->num_contours;
    Contour *contour = &contour_set->contours[contour_id];
    
    // Initialize contour
    contour->num_points = num_points;
    contour->center = contour_set->crack_tip;
    contour->radius = radius;
    contour->contour_type = CONTOUR_CIRCULAR;
    strncpy(contour->name, name, sizeof(contour->name) - 1);
    contour->name[sizeof(contour->name) - 1] = '\0';
    
    // Allocate memory for contour points
    contour->points = (ContourPoint*)malloc(num_points * sizeof(ContourPoint));
    if (!contour->points) {
        fprintf(stderr, "Error: Memory allocation failed for contour points\n");
        return -1;
    }
    
    // Generate circular contour points (counter-clockwise)
    for (int i = 0; i < num_points; i++) {
        double theta = 2.0 * M_PI * i / num_points;
        
        contour->points[i].point.x = contour->center.x + radius * cos(theta);
        contour->points[i].point.y = contour->center.y + radius * sin(theta);
        
        // Outward normal vector (pointing away from center)
        contour->points[i].normal.x = cos(theta);
        contour->points[i].normal.y = sin(theta);
        
        // Differential arc length
        contour->points[i].ds = 2.0 * M_PI * radius / num_points;
        
        // Integration weight (for numerical integration)
        contour->points[i].weight = 1.0;
    }
    
    contour_set->num_contours++;
    
    printf("Added circular contour '%s': radius=%.3f, %d points\n", 
           name, radius, num_points);
    
    return 0;
}

int add_rectangular_contour(ContourSet *contour_set, double width, double height,
                           int num_points_per_side, const char *name) {
    if (!contour_set || width <= 0 || height <= 0 || num_points_per_side < 2) {
        fprintf(stderr, "Error: Invalid parameters for rectangular contour\n");
        return -1;
    }
    
    int contour_id = contour_set->num_contours;
    Contour *contour = &contour_set->contours[contour_id];
    
    int total_points = 4 * num_points_per_side;
    
    // Initialize contour
    contour->num_points = total_points;
    contour->center = contour_set->crack_tip;
    contour->radius = sqrt(width * width + height * height) / 2.0; // Diagonal half-length
    contour->contour_type = CONTOUR_RECTANGULAR;
    strncpy(contour->name, name, sizeof(contour->name) - 1);
    contour->name[sizeof(contour->name) - 1] = '\0';
    
    // Allocate memory for contour points
    contour->points = (ContourPoint*)malloc(total_points * sizeof(ContourPoint));
    if (!contour->points) {
        fprintf(stderr, "Error: Memory allocation failed for contour points\n");
        return -1;
    }
    
    // Rectangle corners relative to center
    double half_width = width / 2.0;
    double half_height = height / 2.0;
    
    int point_idx = 0;
    
    // Bottom side (left to right)
    for (int i = 0; i < num_points_per_side; i++) {
        double t = (double)i / (num_points_per_side - 1);
        contour->points[point_idx].point.x = contour->center.x + (-half_width + t * width);
        contour->points[point_idx].point.y = contour->center.y - half_height;
        contour->points[point_idx].normal.x = 0.0;
        contour->points[point_idx].normal.y = -1.0; // Outward normal (downward)
        contour->points[point_idx].ds = width / num_points_per_side;
        contour->points[point_idx].weight = 1.0;
        point_idx++;
    }
    
    // Right side (bottom to top)
    for (int i = 0; i < num_points_per_side; i++) {
        double t = (double)i / (num_points_per_side - 1);
        contour->points[point_idx].point.x = contour->center.x + half_width;
        contour->points[point_idx].point.y = contour->center.y + (-half_height + t * height);
        contour->points[point_idx].normal.x = 1.0; // Outward normal (rightward)
        contour->points[point_idx].normal.y = 0.0;
        contour->points[point_idx].ds = height / num_points_per_side;
        contour->points[point_idx].weight = 1.0;
        point_idx++;
    }
    
    // Top side (right to left)
    for (int i = 0; i < num_points_per_side; i++) {
        double t = (double)i / (num_points_per_side - 1);
        contour->points[point_idx].point.x = contour->center.x + (half_width - t * width);
        contour->points[point_idx].point.y = contour->center.y + half_height;
        contour->points[point_idx].normal.x = 0.0;
        contour->points[point_idx].normal.y = 1.0; // Outward normal (upward)
        contour->points[point_idx].ds = width / num_points_per_side;
        contour->points[point_idx].weight = 1.0;
        point_idx++;
    }
    
    // Left side (top to bottom)
    for (int i = 0; i < num_points_per_side; i++) {
        double t = (double)i / (num_points_per_side - 1);
        contour->points[point_idx].point.x = contour->center.x - half_width;
        contour->points[point_idx].point.y = contour->center.y + (half_height - t * height);
        contour->points[point_idx].normal.x = -1.0; // Outward normal (leftward)
        contour->points[point_idx].normal.y = 0.0;
        contour->points[point_idx].ds = height / num_points_per_side;
        contour->points[point_idx].weight = 1.0;
        point_idx++;
    }
    
    contour_set->num_contours++;
    
    printf("Added rectangular contour '%s': %.3fx%.3f, %d points\n", 
           name, width, height, total_points);
    
    return 0;
}

int add_custom_contour(ContourSet *contour_set, Point2D *points, 
                      int num_points, const char *name) {
    if (!contour_set || !points || num_points < 3) {
        fprintf(stderr, "Error: Invalid parameters for custom contour\n");
        return -1;
    }
    
    int contour_id = contour_set->num_contours;
    Contour *contour = &contour_set->contours[contour_id];
    
    // Initialize contour
    contour->num_points = num_points;
    contour->center = contour_set->crack_tip;
    contour->contour_type = CONTOUR_CUSTOM;
    strncpy(contour->name, name, sizeof(contour->name) - 1);
    contour->name[sizeof(contour->name) - 1] = '\0';
    
    // Allocate memory for contour points
    contour->points = (ContourPoint*)malloc(num_points * sizeof(ContourPoint));
    if (!contour->points) {
        fprintf(stderr, "Error: Memory allocation failed for contour points\n");
        return -1;
    }
    
    // Copy points and calculate properties
    for (int i = 0; i < num_points; i++) {
        contour->points[i].point = points[i];
        contour->points[i].weight = 1.0;
    }
    
    // Calculate normals and arc lengths
    calculate_contour_normals(contour);
    calculate_contour_arc_lengths(contour);
    
    // Calculate average radius
    double total_radius = 0.0;
    for (int i = 0; i < num_points; i++) {
        total_radius += point_distance(contour->points[i].point, contour->center);
    }
    contour->radius = total_radius / num_points;
    
    contour_set->num_contours++;
    
    printf("Added custom contour '%s': %d points, avg radius=%.3f\n", 
           name, num_points, contour->radius);
    
    return 0;
}

void calculate_contour_normals(Contour *contour) {
    if (!contour || !contour->points) return;
    
    for (int i = 0; i < contour->num_points; i++) {
        int next = (i + 1) % contour->num_points;
        int prev = (i - 1 + contour->num_points) % contour->num_points;
        
        // Calculate tangent vector (average of adjacent segments)
        double tx1 = contour->points[next].point.x - contour->points[i].point.x;
        double ty1 = contour->points[next].point.y - contour->points[i].point.y;
        double tx2 = contour->points[i].point.x - contour->points[prev].point.x;
        double ty2 = contour->points[i].point.y - contour->points[prev].point.y;
        
        double tx = (tx1 + tx2) / 2.0;
        double ty = (ty1 + ty2) / 2.0;
        
        // Normalize tangent
        double t_mag = sqrt(tx * tx + ty * ty);
        if (t_mag > 1e-10) {
            tx /= t_mag;
            ty /= t_mag;
        }
        
        // Outward normal (rotate tangent 90 degrees clockwise)
        contour->points[i].normal.x = ty;
        contour->points[i].normal.y = -tx;
    }
}

void calculate_contour_arc_lengths(Contour *contour) {
    if (!contour || !contour->points) return;
    
    for (int i = 0; i < contour->num_points; i++) {
        int next = (i + 1) % contour->num_points;
        
        double dx = contour->points[next].point.x - contour->points[i].point.x;
        double dy = contour->points[next].point.y - contour->points[i].point.y;
        
        contour->points[i].ds = sqrt(dx * dx + dy * dy);
    }
}

int verify_contour_orientation(const Contour *contour) {
    if (!contour || !contour->points || contour->num_points < 3) return -1;
    
    // Calculate signed area using shoelace formula
    double signed_area = 0.0;
    
    for (int i = 0; i < contour->num_points; i++) {
        int next = (i + 1) % contour->num_points;
        signed_area += (contour->points[next].point.x - contour->points[i].point.x) *
                       (contour->points[next].point.y + contour->points[i].point.y);
    }
    
    // Counter-clockwise if signed area is negative
    return (signed_area < 0) ? 1 : 0;
}

void reverse_contour_orientation(Contour *contour) {
    if (!contour || !contour->points) return;
    
    // Reverse point order
    for (int i = 0; i < contour->num_points / 2; i++) {
        int j = contour->num_points - 1 - i;
        
        ContourPoint temp = contour->points[i];
        contour->points[i] = contour->points[j];
        contour->points[j] = temp;
    }
    
    // Reverse normal directions
    for (int i = 0; i < contour->num_points; i++) {
        contour->points[i].normal.x = -contour->points[i].normal.x;
        contour->points[i].normal.y = -contour->points[i].normal.y;
    }
    
    printf("Contour orientation reversed\n");
}

int point_inside_contour(const Contour *contour, Point2D point) {
    if (!contour || !contour->points) return 0;
    
    // Ray casting algorithm
    int crossings = 0;
    
    for (int i = 0; i < contour->num_points; i++) {
        int next = (i + 1) % contour->num_points;
        
        Point2D p1 = contour->points[i].point;
        Point2D p2 = contour->points[next].point;
        
        // Check if ray crosses edge
        if (((p1.y > point.y) != (p2.y > point.y)) &&
            (point.x < (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y) + p1.x)) {
            crossings++;
        }
    }
    
    return (crossings % 2 == 1);
}

int find_closest_contour_point(const Contour *contour, Point2D point) {
    if (!contour || !contour->points) return -1;
    
    int closest_idx = 0;
    double min_distance = point_distance(contour->points[0].point, point);
    
    for (int i = 1; i < contour->num_points; i++) {
        double distance = point_distance(contour->points[i].point, point);
        if (distance < min_distance) {
            min_distance = distance;
            closest_idx = i;
        }
    }
    
    return closest_idx;
}

int refine_contour(Contour *contour, int refinement_factor) {
    if (!contour || !contour->points || refinement_factor < 2) return -1;
    
    int new_num_points = contour->num_points * refinement_factor;
    ContourPoint *new_points = (ContourPoint*)malloc(new_num_points * sizeof(ContourPoint));
    
    if (!new_points) {
        fprintf(stderr, "Error: Memory allocation failed for refined contour\n");
        return -1;
    }
    
    // Interpolate new points
    for (int i = 0; i < contour->num_points; i++) {
        int next = (i + 1) % contour->num_points;
        
        for (int j = 0; j < refinement_factor; j++) {
            int new_idx = i * refinement_factor + j;
            double t = (double)j / refinement_factor;
            
            // Linear interpolation
            new_points[new_idx].point.x = (1.0 - t) * contour->points[i].point.x + 
                                         t * contour->points[next].point.x;
            new_points[new_idx].point.y = (1.0 - t) * contour->points[i].point.y + 
                                         t * contour->points[next].point.y;
            
            new_points[new_idx].normal.x = (1.0 - t) * contour->points[i].normal.x + 
                                          t * contour->points[next].normal.x;
            new_points[new_idx].normal.y = (1.0 - t) * contour->points[i].normal.y + 
                                          t * contour->points[next].normal.y;
            
            new_points[new_idx].ds = contour->points[i].ds / refinement_factor;
            new_points[new_idx].weight = 1.0;
        }
    }
    
    // Replace old points
    free(contour->points);
    contour->points = new_points;
    contour->num_points = new_num_points;
    
    printf("Contour refined: %d points -> %d points\n", 
           contour->num_points / refinement_factor, contour->num_points);
    
    return 0;
}

void print_contour_info(const Contour *contour) {
    if (!contour) {
        printf("Error: Null contour pointer\n");
        return;
    }
    
    printf("\n=== CONTOUR INFORMATION ===\n");
    printf("Name: %s\n", contour->name);
    printf("Type: %s\n", 
           (contour->contour_type == CONTOUR_CIRCULAR) ? "Circular" :
           (contour->contour_type == CONTOUR_RECTANGULAR) ? "Rectangular" : "Custom");
    printf("Center: (%.6f, %.6f)\n", contour->center.x, contour->center.y);
    printf("Radius: %.6f\n", contour->radius);
    printf("Number of points: %d\n", contour->num_points);
    
    // Calculate total contour length
    double total_length = 0.0;
    for (int i = 0; i < contour->num_points; i++) {
        total_length += contour->points[i].ds;
    }
    printf("Total length: %.6f\n", total_length);
    
    printf("Orientation: %s\n", 
           verify_contour_orientation(contour) ? "Counter-clockwise" : "Clockwise");
    printf("===========================\n\n");
}

void print_contour_set_info(const ContourSet *contour_set) {
    if (!contour_set) {
        printf("Error: Null contour set pointer\n");
        return;
    }
    
    printf("\n=== CONTOUR SET INFORMATION ===\n");
    printf("Crack tip: (%.6f, %.6f)\n", contour_set->crack_tip.x, contour_set->crack_tip.y);
    printf("Number of contours: %d\n", contour_set->num_contours);
    
    for (int i = 0; i < contour_set->num_contours; i++) {
        printf("\nContour %d:\n", i);
        print_contour_info(&contour_set->contours[i]);
    }
    printf("===============================\n\n");
}

void free_contour(Contour *contour) {
    if (!contour) return;
    
    if (contour->points) {
        free(contour->points);
        contour->points = NULL;
    }
}

void free_contour_set(ContourSet *contour_set) {
    if (!contour_set) return;
    
    if (contour_set->contours) {
        for (int i = 0; i < contour_set->num_contours; i++) {
            free_contour(&contour_set->contours[i]);
        }
        free(contour_set->contours);
    }
    
    free(contour_set);
}

int validate_contour(const Contour *contour) {
    if (!contour || !contour->points) return -1;
    
    int errors = 0;
    
    // Check for valid points
    for (int i = 0; i < contour->num_points; i++) {
        if (!isfinite(contour->points[i].point.x) || !isfinite(contour->points[i].point.y) ||
            !isfinite(contour->points[i].normal.x) || !isfinite(contour->points[i].normal.y) ||
            !isfinite(contour->points[i].ds) || contour->points[i].ds <= 0) {
            errors++;
        }
    }
    
    if (errors > 0) {
        printf("Warning: Found %d invalid contour points\n", errors);
    }
    
    return errors;
} 