import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class CellularAutomata:
    def __init__(self, N, p1, p2):
        self.N = N
        self.p1 = p1
        self.p2 = p2
        self.sweeps_per_update = 1
        # Init a grid of size NxN with 0=R, 1=G, 2=B with p=1/3 each.
        self.grid = np.random.choice([0, 1, 2], p=[1/3, 1/3, 1/3], size=(N, N))

        self.animation = True

    def choose_nn(self, i, j):
        choice = np.random.randint(low=0, high=3)
        if choice == 0:
            return self.grid[(i+1+self.N)%self.N][j]
        if choice == 1:
            return self.grid[(i-1+self.N)%self.N][j]
        if choice == 2:
            return self.grid[i][(j+1+self.N)%self.N]
        if choice == 3:
            return self.grid[i][(i-1+self.N)%self.N]

    def update(self, k):
        for z in range(self.sweeps_per_update):
            rand_choices = np.random.randint(low=0, high=self.N, size=(self.N, self.N, 2))
            update_probabilities = np.random.uniform(size=(self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    # print(rand_choices[i][j][0], rand_choices[i][j][1])
                    print(update_probabilities[i][j])
                    val = self.choose_nn(rand_choices[i][j][0], rand_choices[i][j][1])
                    # since the update is cyclic RG -> GB -> BR = p1 and the reverse p2
                    # if val > the selected cell use p1 else use p2.
                    if val > self.grid[rand_choices[i][j][0]][rand_choices[i][j][1]] and update_probabilities[i][j] <= self.p1:
                        self.grid[rand_choices[i][j][0]][rand_choices[i][j][1]] = val
                    elif val < self.grid[rand_choices[i][j][0]][rand_choices[i][j][1]] and update_probabilities[i][j] <= self.p2:
                        self.grid[rand_choices[i][j][0]][rand_choices[i][j][1]] = val
                    else:
                        continue
        if self.animation:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower')# , vmim=-1, vmax=2)
            # plt.colorbar()


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
