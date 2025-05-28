#define _USE_MATH_DEFINES
#include "../include/mesh.h"
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Mesh* create_rectangular_mesh(int nx, int ny, double width, double height) {
    if (nx < 2 || ny < 2 || width <= 0 || height <= 0) {
        fprintf(stderr, "Error: Invalid mesh parameters\n");
        return NULL;
    }
    
    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    if (!mesh) {
        fprintf(stderr, "Error: Memory allocation failed for mesh\n");
        return NULL;
    }
    
    mesh->nx = nx;
    mesh->ny = ny;
    mesh->num_nodes = nx * ny;
    mesh->num_elements = (nx - 1) * (ny - 1);
    mesh->width = width;
    mesh->height = height;
    
    // Allocate memory for nodes
    mesh->nodes = (Node*)malloc(mesh->num_nodes * sizeof(Node));
    if (!mesh->nodes) {
        fprintf(stderr, "Error: Memory allocation failed for nodes\n");
        free(mesh);
        return NULL;
    }
    
    // Allocate memory for elements
    mesh->elements = (Element*)malloc(mesh->num_elements * sizeof(Element));
    if (!mesh->elements) {
        fprintf(stderr, "Error: Memory allocation failed for elements\n");
        free(mesh->nodes);
        free(mesh);
        return NULL;
    }
    
    // Generate nodes
    double dx = width / (nx - 1);
    double dy = height / (ny - 1);
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int node_id = j * nx + i;
            mesh->nodes[node_id].id = node_id;
            mesh->nodes[node_id].coord.x = i * dx;
            mesh->nodes[node_id].coord.y = j * dy;
        }
    }
    
    // Generate elements (quadrilaterals)
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int elem_id = j * (nx - 1) + i;
            
            // Node connectivity (counter-clockwise)
            mesh->elements[elem_id].nodes[0] = j * nx + i;           // bottom-left
            mesh->elements[elem_id].nodes[1] = j * nx + (i + 1);     // bottom-right
            mesh->elements[elem_id].nodes[2] = (j + 1) * nx + (i + 1); // top-right
            mesh->elements[elem_id].nodes[3] = (j + 1) * nx + i;     // top-left
            
            // Calculate element area
            mesh->elements[elem_id].area = dx * dy;
        }
    }
    
    // Initialize crack (no crack initially)
    mesh->crack.start.x = 0.0;
    mesh->crack.start.y = 0.0;
    mesh->crack.end.x = 0.0;
    mesh->crack.end.y = 0.0;
    mesh->crack.tip.x = 0.0;
    mesh->crack.tip.y = 0.0;
    mesh->crack.length = 0.0;
    mesh->crack.angle = 0.0;
    
    return mesh;
}

void define_crack(Mesh *mesh, double start_x, double start_y, 
                  double end_x, double end_y) {
    if (!mesh) {
        fprintf(stderr, "Error: Null mesh pointer\n");
        return;
    }
    
    mesh->crack.start.x = start_x;
    mesh->crack.start.y = start_y;
    mesh->crack.end.x = end_x;
    mesh->crack.end.y = end_y;
    
    // Calculate crack length
    double dx = end_x - start_x;
    double dy = end_y - start_y;
    mesh->crack.length = sqrt(dx * dx + dy * dy);
    
    // Calculate crack angle (in radians)
    mesh->crack.angle = atan2(dy, dx);
    
    // Set crack tip (assuming crack grows from start to end)
    mesh->crack.tip.x = end_x;
    mesh->crack.tip.y = end_y;
    
    printf("Crack defined: start(%.3f, %.3f) -> end(%.3f, %.3f)\n",
           start_x, start_y, end_x, end_y);
    printf("Crack length: %.3f, angle: %.3f rad (%.1f deg)\n",
           mesh->crack.length, mesh->crack.angle, mesh->crack.angle * 180.0 / M_PI);
}

Mesh* load_mesh_from_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open mesh file %s\n", filename);
        return NULL;
    }
    
    int nx, ny;
    double width, height;
    
    // Read mesh parameters
    if (fscanf(file, "%d %d", &nx, &ny) != 2) {
        fprintf(stderr, "Error: Cannot read mesh dimensions\n");
        fclose(file);
        return NULL;
    }
    
    if (fscanf(file, "%lf %lf", &width, &height) != 2) {
        fprintf(stderr, "Error: Cannot read mesh size\n");
        fclose(file);
        return NULL;
    }
    
    // Create mesh
    Mesh *mesh = create_rectangular_mesh(nx, ny, width, height);
    if (!mesh) {
        fclose(file);
        return NULL;
    }
    
    // Read crack definition if present
    double crack_x1, crack_y1, crack_x2, crack_y2;
    if (fscanf(file, "%lf %lf %lf %lf", &crack_x1, &crack_y1, &crack_x2, &crack_y2) == 4) {
        define_crack(mesh, crack_x1, crack_y1, crack_x2, crack_y2);
    }
    
    fclose(file);
    return mesh;
}

void save_mesh_to_file(const Mesh *mesh, const char *filename) {
    if (!mesh || !filename) {
        fprintf(stderr, "Error: Invalid parameters for save_mesh_to_file\n");
        return;
    }
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Cannot create mesh file %s\n", filename);
        return;
    }
    
    // Write mesh parameters
    fprintf(file, "# Mesh dimensions (nx ny)\n");
    fprintf(file, "%d %d\n", mesh->nx, mesh->ny);
    
    fprintf(file, "# Domain size (width height)\n");
    fprintf(file, "%.6f %.6f\n", mesh->width, mesh->height);
    
    // Write node coordinates
    fprintf(file, "# Node coordinates (x y)\n");
    for (int i = 0; i < mesh->num_nodes; i++) {
        fprintf(file, "%.6f %.6f\n", 
                mesh->nodes[i].coord.x, mesh->nodes[i].coord.y);
    }
    
    // Write crack definition
    fprintf(file, "# Crack definition (start_x start_y end_x end_y)\n");
    fprintf(file, "%.6f %.6f %.6f %.6f\n",
            mesh->crack.start.x, mesh->crack.start.y,
            mesh->crack.end.x, mesh->crack.end.y);
    
    fclose(file);
    printf("Mesh saved to %s\n", filename);
}

int find_element_containing_point(const Mesh *mesh, Point2D point) {
    if (!mesh) return -1;
    
    double dx = mesh->width / (mesh->nx - 1);
    double dy = mesh->height / (mesh->ny - 1);
    
    // Find element indices
    int i = (int)(point.x / dx);
    int j = (int)(point.y / dy);
    
    // Check bounds
    if (i < 0 || i >= mesh->nx - 1 || j < 0 || j >= mesh->ny - 1) {
        return -1;
    }
    
    return j * (mesh->nx - 1) + i;
}

Point2D get_node_coordinates(const Mesh *mesh, int node_id) {
    Point2D invalid_point = {-1.0, -1.0};
    
    if (!mesh || node_id < 0 || node_id >= mesh->num_nodes) {
        return invalid_point;
    }
    
    return mesh->nodes[node_id].coord;
}

double point_distance(Point2D p1, Point2D p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return sqrt(dx * dx + dy * dy);
}

void free_mesh(Mesh *mesh) {
    if (!mesh) return;
    
    if (mesh->nodes) {
        free(mesh->nodes);
    }
    
    if (mesh->elements) {
        free(mesh->elements);
    }
    
    free(mesh);
}

void print_mesh_info(const Mesh *mesh) {
    if (!mesh) {
        printf("Error: Null mesh pointer\n");
        return;
    }
    
    printf("\n=== MESH INFORMATION ===\n");
    printf("Dimensions: %d x %d nodes\n", mesh->nx, mesh->ny);
    printf("Total nodes: %d\n", mesh->num_nodes);
    printf("Total elements: %d\n", mesh->num_elements);
    printf("Domain size: %.3f x %.3f\n", mesh->width, mesh->height);
    
    if (mesh->crack.length > 0) {
        printf("\n--- CRACK INFORMATION ---\n");
        printf("Start: (%.3f, %.3f)\n", mesh->crack.start.x, mesh->crack.start.y);
        printf("End: (%.3f, %.3f)\n", mesh->crack.end.x, mesh->crack.end.y);
        printf("Tip: (%.3f, %.3f)\n", mesh->crack.tip.x, mesh->crack.tip.y);
        printf("Length: %.3f\n", mesh->crack.length);
        printf("Angle: %.3f rad (%.1f deg)\n", 
               mesh->crack.angle, mesh->crack.angle * 180.0 / M_PI);
    } else {
        printf("No crack defined\n");
    }
    printf("========================\n\n");
} 