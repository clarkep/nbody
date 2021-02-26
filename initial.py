import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# Define the parameters for the problem
N = 10     # Number of bodies
X_MIN = -100    # Coordinate boundaries for the initial positions of the bodies
X_MAX = 100
Y_MIN = -100
Y_MAX = 100
Z_MIN = 0
Z_MAX = 0


# For now, this just randomly places the bodies in space.
def main():
    bodies = np.zeros((N, 7))
    X, Y, Z, VX, VY, VZ, MASS = range(7)
    bodies[:, X] = np.random.uniform(X_MIN, X_MAX, size=N)
    bodies[:, Y] = np.random.uniform(Y_MIN, Y_MAX, size=N)
    bodies[:, Z] = np.random.uniform(Z_MIN, Z_MAX, size=N)

    # Give small initial velocities to break symmetry
    bodies[:, VX] = np.random.uniform(-1, 1, size=N)
    bodies[:, VY] = np.random.uniform(-1, 1, size=N)
    # bodies[:, VZ] = np.random.uniform(-1, 1, size=N)

    bodies[:, MASS] = np.random.uniform(0.1, 20, size=N)

    bodies = pd.DataFrame(bodies)
    bodies.columns = 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'
    bodies.to_csv('products/init.csv', index=False)


if __name__ == '__main__':
    main()

