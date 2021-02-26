#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 10   // number of bodies. Make sure this agrees with the input file!
#define INIT_FILE "products/init.csv"
#define OUTPUT_FILE "products/output.csv"
#define TIME_STEP 0.1
#define DURATION 1000.0
#define OUTPUT_EVERY 10     // time steps between writing to the output file
#define G 1.0   // gravitational constant
#define EPS 2.    // Smoothing constant to avoid singularities in the forces


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
        printf("%f,%d,%f,%f,%f,%f,%f,%f,%f\n",
               t, i, bodies[i].x, bodies[i].y, bodies[i].z,
                     bodies[i].vx, bodies[i].vy, bodies[i].vz, bodies[i].mass);
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

        // TODO: Could we use Runge-Kutta here instead of Euler? Could be tricky
        //  to get the "value of the function" for t + 0.5dt, etc. We would have
        //  to simulate the bodies for several different possible updates, and
        //  average the forces, or something. Weird. Ask Sasha.
        // See: https://www.physicsforums.com/threads/runge-kutta-for-gravitational-n-body-simulation-prediction-of-acceleration.599677/
        // Also, which will give us more error: the Barnes-Hut approximation, or Euler approximation?

        bodies[i].vx += accel_x * TIME_STEP;
        bodies[i].vy += accel_y * TIME_STEP;
        bodies[i].vz += accel_z * TIME_STEP;
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












