#include <stdlib.h>
#include <stdio.h>
#include <math.h>


typedef struct {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double mass;
} Body;

// Node for the Barnes-Hut tree. All coordinates are absolute (not relative
// to the node's position).
struct Node {
    struct Node **children;    // Array of pointers to this node's children

    double side_length;
    // Coordinates of the center of the node's region
    double x;
    double y;

    // Coordinates of the center of mass of the bodies in the node. If there is
    // just one body, then this is just the position of the body.
    double center_of_mass_x;
    double center_of_mass_y;

    double total_mass;
};
typedef struct Node Node;

// Reads the input CSV file of body positions and velocities, and returns an
// array of Bodies
Body* read_init(char *filename) {

    // Allocate space on the heap for the bodies array
    Body* bodies = (Body*) malloc(N * sizeof(Body));

    FILE *init = fopen(filename, "r");
    // Throw away the header line
    fscanf(init, "%*s");

    int i = 0;
    float x, y, z, vx, vy, vz, mass;
    while (fscanf(init, "%f,%f,%f,%f,%f,%f,%f",
                        &x, &y, &z, &vx, &vy, &vz, &mass) != EOF) {
        bodies[i].x = x;
        bodies[i].y = y;
        bodies[i].z = z;
        bodies[i].vx = vx;
        bodies[i].vy = vy;
        bodies[i].vz = vz;
        bodies[i].mass = mass;
        i += 1;
    }

    return bodies;
}

// Print out the current state of the system (the time, and positions of
// all the bodies and their masses)
void write_state(Body* bodies, double t) {
    for (int i = 0; i < N; i++) {
        printf("%f,%d,%f,%f,%f,%f,%f,%f,%f\n",
               t, i, bodies[i].x, bodies[i].y, bodies[i].z,
                     bodies[i].vx, bodies[i].vy, bodies[i].vz, bodies[i].mass);
    }
}

// Update the positions and velocities of each body based on the gravitational
// forces on it.
void naive_time_step(Body* bodies) {
    for (int i = 0; i < N; i++) {

        // Update the position based on the velocity
        bodies[i].x += bodies[i].vx * TIME_STEP;
        bodies[i].y += bodies[i].vy * TIME_STEP;
        bodies[i].z += bodies[i].vz * TIME_STEP;

        // Update the velocity based on the acceleration from the gravitational forces
        double m1 = bodies[i].mass;
        double force_x = 0;
        double force_y = 0;
        double force_z = 0;   // net force in each direction
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double m2 = bodies[j].mass;
            double x1 = bodies[i].x;
            double x2 = bodies[j].x;
            double y1 = bodies[i].y;
            double y2 = bodies[j].y;
            double z1 = bodies[i].z;
            double z2 = bodies[j].z;

            // See the documentation for the formulas here
            double r2 = pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2) + EPS;
            double r3 = pow(r2, 1.5);

            double force = G * m1 * m2 / r3;
            force_x += force * (x2 - x1);
            force_y += force * (y2 - y1);
            force_z += force * (z2 - z1);
        }

        double accel_x = force_x / m1;
        double accel_y = force_y / m1;
        double accel_z = force_z / m1;

        bodies[i].vx += accel_x * TIME_STEP;
        bodies[i].vy += accel_y * TIME_STEP;
        bodies[i].vz += accel_z * TIME_STEP;
    }
}


// If n already has its center position and side length fields, then fill out the center of
// mass and total mass. This iterates through bodies to find the center of mass
// and total mass of the bodies inside the node. Note: If there are no bodies in
// the node's region, then the center of mass and total mass will all be set to 0.
void set_mass_info(Node *n, Body *bodies) {
    double cm_x = 0.0;
    double cm_y = 0.0;
    double total_mass = 0.0;

    double max_x = n->x + n->side_length / 2.0;
    double min_x = n->x - n->side_length / 2.0;
    double max_y = n->y + n->side_length / 2.0;
    double min_y = n->y - n->side_length / 2.0;

    for (int i=0; i<N; i++) {
        if (bodies[i].x > min_x
            && bodies[i].x < max_x
            && bodies[i].y > min_y
            && bodies[i].y < max_y) {

                total_mass += bodies[i].mass;
                cm_x += bodies[i].mass * bodies[i].x;
                cm_y += bodies[i].mass * bodies[i].y;
        }
    }

    n->center_of_mass_x = cm_x / total_mass;
    n->center_of_mass_y = cm_y / total_mass;
    n->total_mass = total_mass;
}

// Alternate design idea: each node stores an array of the bodies it contains.
// Do this if the O(N) checking for bodies inside the node takes too long.

// Returns the number of bodies that are inside the node's region.
int bodies_in_node(Node *n, Body *bodies) {
    int counter = 0;
    double max_x = n->x + n->side_length / 2.0;
    double min_x = n->x - n->side_length / 2.0;
    double max_y = n->y + n->side_length / 2.0;
    double min_y = n->y - n->side_length / 2.0;

    for (int i=0; i<N; i++) {
        if (bodies[i].x > min_x
            && bodies[i].x < max_x
            && bodies[i].y > min_y
            && bodies[i].y < max_y) {
                counter += 1;
        }
    }

    return counter;
}


void build_subtree(Node *n, Body *bodies) {
    if (bodies_in_node(n, bodies) > 1) {
        n->children = (Node **) malloc(4 * sizeof(Node *));
        for (int i=0; i<4; i++) {
            Node *child = malloc(sizeof(Node));
            child->children = NULL;
            child->side_length = n->side_length / 2.0;

            // Children as arranged like this:
            // 0 1
            // 2 3
            if (i == 0) {
                child->x = n->x - child->side_length / 2.0;
                child->y = n->y + child->side_length / 2.0;
            } else if (i == 1) {
                child->x = n->x + child->side_length / 2.0;
                child->y = n->y + child->side_length / 2.0;
            } else if (i == 2) {
                child->x = n->x - child->side_length / 2.0;
                child->y = n->y - child->side_length / 2.0;
            } else if (i == 3) {
                child->x = n->x + child->side_length / 2.0;
                child->y = n->y - child->side_length / 2.0;
            }

            // fprintf(stderr, "(%f, %f, %f),\n", child->x, child->y, child->side_length);

            set_mass_info(child, bodies);
            build_subtree(child, bodies);
            n->children[i] = child;
        }
    }
}


// Figure out how big the side length should be of the root node, which is 2x the
// max x, y, or z value amongst all the bodies.
double max_side_length(Body *bodies) {
    double max = 0.0;
    for (int i=0; i<N; i++) {
        if (bodies[i].x > max) {
            max = bodies[i].x;
        }
        if (bodies[i].y > max) {
            max = bodies[i].y;
        }
        if (bodies[i].z > max) {
            max = bodies[i].z;
        }
    }
    return 2*max;
}

// Given an array of Bodies, construct a Barnes-Hut tree and return the pointer
// to the root node.
Node *barnes_hut_tree(Body *bodies) {
    Node *n = malloc(sizeof(Node));
    n->side_length = 200.0; // max_side_length(bodies);
    // The center of the root node is the origin
    n->x = 0.0;
    n->y = 0.0;
    set_mass_info(n, bodies);
    build_subtree(n, bodies);
    return n;
}

// Returns the net gravitational force on the given body (as an array of {x, y, z})
void net_force(Body body, Node *tree, double *force) {

    if (tree == NULL) return;

    double distance = sqrt(pow(body.x-tree->center_of_mass_x, 2.)
                           + pow(body.y-tree->center_of_mass_y, 2.)
                           // + pow(body.z-tree->center_of_mass_z, 2.)
                           + EPS);     // Epsilon smoothing to avoid singularities

    if (tree->side_length / distance < THETA || tree->children == NULL) {
        // Base case: s/d < theta
        // I'm doing the thing where you divide by r^3 instead of r^2 so that
        // you multiply by \vec{r} instead of \hat{r}
        double force_mag = G * body.mass * tree->total_mass / pow(distance, 3.0);
        force[0] += force_mag * (tree->center_of_mass_x - body.x);
        force[1] += force_mag * (tree->center_of_mass_y - body.y);
        // force[2] += force_mag * (tree->center_of_mass_z - body.z);
    } else {
        // Recursive case: go down into the node's children (sub-squares) and
        //  find the force from those.
        for (int i=0; i<4; i++) {
            net_force(body, tree->children[i], force);
        }
    }
}

void barnes_hut_step(Body *bodies) {
    Node *tree = barnes_hut_tree(bodies);
    for (int i=0; i<N; i++) {
        double force[3] = {0.0, 0.0, 0.0};
        net_force(bodies[i], tree, force);
        double accel_x = force[0] / bodies[i].mass;
        double accel_y = force[1] / bodies[i].mass;
        double accel_z = force[2] / bodies[i].mass;

        // update the velocity based on the acceleration
        bodies[i].vx += accel_x * TIME_STEP;
        bodies[i].vy += accel_y * TIME_STEP;
        bodies[i].vz += accel_z * TIME_STEP;

        // update the position based on the velocity
        bodies[i].x += bodies[i].vx * TIME_STEP;
        bodies[i].y += bodies[i].vy * TIME_STEP;
        bodies[i].z += bodies[i].vz * TIME_STEP;
    }
}


// Naive O(n^2) solver
int main(int argc, char *argv[]) {
    Body *bodies = read_init(argv[1]);
    int max_time_steps = (int) DURATION / TIME_STEP;
    for (int i=0; i<max_time_steps; i++) {
        barnes_hut_step(bodies);
        // naive_time_step(bodies);
        if (i % OUTPUT_EVERY == 0) {
            double t = i * TIME_STEP;
            write_state(bodies, t);
        }
        fprintf(stderr, "\rTime step: %d/%d", i+1, max_time_steps);
        fflush(stderr);
    }
}












