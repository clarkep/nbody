import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from graphics import *
import time

SCALE = 1
WIDTH = 800
HEIGHT = 600
N = 10


def clear(win):
    for item in win.items[:]:
        item.undraw()


def animate():
    df = pd.read_csv('products/output.csv', header=None)
    df.columns = 't', 'id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'

    win = GraphWin('Pendulum', WIDTH, HEIGHT, autoflush=False) # give title and dimensions
    for i in range(0, df.shape[0], 1):

        if i % N == 0:
            clear(win)

        x = df.loc[i]['x'] * SCALE + WIDTH/2
        y = -df.loc[i]['y'] * SCALE + HEIGHT/2
        mass = df.loc[i]['mass']
        id = df.loc[i]['id']
        # origin = Point(WIDTH/2, HEIGHT/2)
        location = Point(x, y)
        body = Circle(Point(x, y), mass/2)
        body.setFill('black')
        body.draw(win)

        if id == N-1:
            win.update()


    win.close()


if __name__ == '__main__':
    animate()
