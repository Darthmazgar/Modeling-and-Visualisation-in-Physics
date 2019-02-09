import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class GameOfLife:
    def __init__(self, N, M, anim=True):
        self.N = N
        self.M = M
        self.sweeps = 1
        self.grid = np.random.choice([0, 1], size=(N, M), p=[.9, .1])
        self.new_grid = self.grid
        if anim:
            self.fig = plt.figure()
        # self.init_kaw_grid()

    def init_kaw_grid(self):
        """
        Initialise the grid with one hald all spin 1 and the other all spin -1.
        Used for initialising Kawasaki tests.
        """
        ones = np.ones((int(self.N / 2), self.M))
        neg_ones = ones * -1
        self.grid = np.concatenate((ones, neg_ones))
        return self.grid

    def update(self, k, anim=True):
        if anim:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='None',
                       cmap='Blues', vmin=0, vmax=1)
        for z in range(self.sweeps):
            for i in range(self.N):
                for j in range(self.M):
                    count = 0
                    state = self.grid[i][j]
                    for x in range(-1, 2):
                        for y in range(-1, 2):
                            if x == 0 and y == 0:
                                continue
                            count += self.grid[(i + x) % self.N]\
                                              [(j + y) % self.M]

                    # count -= state
                    if state == 0 and count == 3:
                        self.new_grid[i][j] = 1
                    elif state == 1 and (count < 2 or count > 3):
                        self.new_grid[i][j] = 0
            self.grid = self.new_grid

        return self.grid

    def population(self, pr=True):
        size = self.N * self.M
        pop = np.sum(self.grid)
        frac = pop / size
        if pr:
            print("Alive population: %.1f%%" % frac*100)
        return frac

    def animate(self):
        anim = FuncAnimation(self.fig, self.update)
        plt.show()


def main(argv):
    if len(argv) != 2:
        print("gol.py N M")
        sys.exit()
    N = int(argv[0])
    M = int(argv[1])
    gol = GameOfLife(N, M)
    gol.animate()

main(sys.argv[1:])
