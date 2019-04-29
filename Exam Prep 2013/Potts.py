import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import exp
import sys

class PottsModel:
    def __init__(self, N, T, kb=1, J=1):
        self.N = N
        self.T = T
        self.kb = kb
        self.J = J
        self.grid = np.random.choice([0,1,2], size=(N,N))
        self.animation = True
        self.sweeps_per_update = 1

    def update(self, k):
        for z in range(self.sweeps_per_update):
            choices = np.random.randint(0, self.N, size=(self.N, self.N, 2))
            update_probabilities = np.random.uniform(size=(self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    E = self.E(choices[i][j][0], choices[i][j][1])
                    # E = self.E(i, j)
                    print("E", E)
                    P = self.P(E)
                    # P = 0.5
                    print("P", P)
                    if update_probabilities[i][j] > P:

                        val = self.grid[choices[i][j][0]][choices[i][j][1]]
                        # val = self.grid[i][j]
                        self.grid[choices[i][j][0]][choices[i][j][1]] = np.random.choice([(val+1+3)%3, (val-1+3)%3])
                        # self.grid[i][j] = np.random.choice([(val+1+3)%3, (val-1+3)%3])
                        # self.grid[i][j] = 2# np.random.choice([(val+1+3)%3, (val-1+3)%3])
                    else:
                        # print("Not updated")
                        pass
        if self.animation:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower' , vmin=0, vmax=2)
            plt.colorbar()

    def P(self, E):
        return exp(-E / (self.kb * self.T))

    def E(self, i, j):
        home = self.grid[i][j]
        count = 0
        if self.grid[(i+1+self.N) % self.N][j] == home:
            count += 1
        if self.grid[(i-1+self.N) % self.N][j] == home:
            count += 1
        if self.grid[i][(j+1+self.N) % self.N] == home:
            count += 1
        if self.grid[i][(j+1+self.N) % self.N] == home:
            count += 1
        return -self.J * count

    def run_animation(self):
        """
        Gives the ability to click on the animation canvas to play and pause.
        """
        self.animation = True
        self.fig = plt.figure()
        anim_running = True

        def onClick(event):
            nonlocal anim_running
            if anim_running:
                anim.event_source.stop()
                anim_running = False
            else:
                anim.event_source.start()
                anim_running = True

        self.fig.canvas.mpl_connect('button_press_event', onClick)
        anim = FuncAnimation(self.fig, self.update, interval=5)
        plt.show()
