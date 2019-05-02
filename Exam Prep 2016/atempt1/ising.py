import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import sys


class Ising:
    def __init__(self, N=50, h=0.01, J=-1, kb=1, T=1):
        self.N = N
        self.h = h
        self.J = J
        self.kb = kb
        self.T = T
        self.grid = np.random.choice([-1, 1], p=[.5, .5], size=(N, N))
        self.animation = True
        self.sweeps_per_update = 1

    def ising_update(self, coords):
        i = coords[0]
        j = coords[1]
        sum_i_j = (self.grid[(i+1+self.N) % self.N][j]
                + self.grid[(i-1+self.N) % self.N][j]
                + self.grid[i][(j+1+self.N) % self.N]
                + self.grid[i][(j-1+self.N) % self.N])

        Eb = -self.J*sum_i_j - self.h*self.grid[i][j]
        Ea = -self.J*sum_i_j + self.h*self.grid[i][j]
        de = Ea - Eb
        # print(de)
        if de <= 0:
            self.grid[i][j] *= -1
            return True

        elif np.exp(de/(self.kb*self.T)) <= np.random.uniform():
            self.grid[i][j] *= -1
            return True
        else:
            return False

    def update(self, k):
        for z in range(self.sweeps_per_update):
            rand_choices = np.random.randint(self.N, size=(self.N, self.N, 2))
            # print(rand_choices)
            for i in range(self.N):
                for j in range(self.N):
                    # self.grid[rand_choices[i][j][0]][rand_choices[i][j][1]] = self.ising_update(rand_choices[i][j])
                    self.ising_update(rand_choices[i][j])

        if self.animation:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower')# ,
                           # vmim=-1, vmax=1)
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
                print("Paused.")
                anim.event_source.stop()
                anim_running = False
            else:
                print("Resume.")
                anim.event_source.start()
                anim_running = True
        self.fig.canvas.mpl_connect('button_press_event', onClick)
        anim = FuncAnimation(self.fig, self.update, interval=25)
        plt.show()
