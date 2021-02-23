#include <stdlib.h>
#include <stdio.h>

#define N 100   // number of bodies
#define INIT_FILE "products/init.csv"
#define OUTPUT_FILE "products/output.csv"
#define TIME_STEP 0.01
#define DURATION 100.0
#define OUTPUT_EVERY 20     // time steps between writing to the output file

typedef struct {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double mass;
} Body;


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
        printf("%f,%d,%f,%f,%f,%f\n",
               t, i, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].mass);
    }
}

// Update the positions and velocities of each body based on the gravitational
// forces on it.
void time_step(Body* bodies) {
    for (int i = 0; i < N; i++) {

        // Update the position based on the velocity
        bodies[i].x += bodies[i].vx * TIME_STEP;
        bodies[i].y += bodies[i].vy * TIME_STEP;
        bodies[i].z += bodies[i].vz * TIME_STEP;

        // Update the velocity based on the acceleration
        bodies[i].vy -= 0.01;
    }
}


// Naive O(n^2) solver
int main(int argc, char *argv[]) {
    Body *bodies = read_init(argv[1]);
    int max_time_steps = DURATION / TIME_STEP;
    for (int i=0; i<max_time_steps; i++) {
        time_step(bodies);
        if (i % OUTPUT_EVERY == 0) {
            double t = i * TIME_STEP;
            write_state(bodies, t);
        }
    }
}
