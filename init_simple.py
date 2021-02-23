import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# Define the parameters for the problem
X_MIN = -100    # Coordinate boundaries for the initial positions of the bodies
X_MAX = 100
Y_MIN = -100
Y_MAX = 100
Z_MIN = 0
Z_MAX = 0


# For now, this just randomly places the bodies in space.
def main():
    # two bodies, with total momentum=0, total ang. mom. in +z direction.  
    bodies = np.array([[-50, 0, 0, 0, -1, 0, 10], 
                       [ 50, 0, 0, 0,  1, 0, 10]])

    bodies = pd.DataFrame(bodies)
    bodies.columns = 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'
    bodies.to_csv('products/init_simple.csv', index=False)


if __name__ == '__main__':
    main()

