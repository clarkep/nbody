
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# Define the parameters for the problem

N = int(sys.argv[1])     # Number of bodies
R = 50  #Plummer radius

X, Y, Z, VX, VY, VZ, MASS = range(7)

# .
def main():
    bodies = np.zeros((N, 7))
    Radius = R/(np.sqrt(np.power(np.random.rand(N),-2/3)-1)) #Plummer model radius distribution
    Theta = np.pi * np.random.rand(N)   #Randomized spherical coordinate angles
    Phi = 2 * np.pi * np.random.rand(N)
    
    bodies[:, X] = Radius * np.sin(Theta) * np.cos (Phi)
    bodies[:, Y] = Radius * np.sin(Theta) * np.sin (Phi)
    bodies[:, Z] = Radius * np.cos(Theta)

    # Rotation
    bodies[:, VX] = bodies[:, Y] / R 
    bodies[:, VY] = -bodies[:, X] / R 
    bodies[:, VZ] = bodies[:, Y] / R 


    bodies[:, MASS] = np.random.uniform(0.1, 20, size=N)

    bodies = pd.DataFrame(bodies)
    bodies.columns = 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'
    bodies.to_csv('products/init.csv', index=False)


if __name__ == '__main__':
    # figure_eight()
    main()
