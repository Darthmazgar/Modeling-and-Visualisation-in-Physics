import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys


class Poisson:
    def __init__(self, int N, float accuracy, float dx=1, float epsilon=1):
        self.N = N
        self.acc = accuracy
        self.dx = dx
        self.epsilon = epsilon
        self.rho_grid = np.zeros((N, N, N))
        self.phi_grid = np.zeros((N, N, N))
        self.e_grid = None


    def add_point_charge(self, int i, int j, int k):
        if i >= self.N or j >= self.N or k >= self.N:
            raise ValueError("Point selected (%d, %d, %d) out of range" % (i, j, k)
                        + " (%d, %d, %d)" % (self.N, self.N, self.N))
        self.rho_grid[i][j][k] = 1


    def update(self, int k):
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    self.phi_grid[i][j][k] = self.jacobi_update(i, j, k)
                    # self.next_phi(i, j, k)

        self.zero_boundaries()

    def zero_boundaries(self):
        self.phi_grid[:, :, 0] = 0
        self.phi_grid[:, :, -1] = 0
        self.phi_grid[:, 0, :] = 0
        self.phi_grid[:, -1, :] = 0
        self.phi_grid[0, :, :] = 0
        self.phi_grid[-1, :, :] = 0

    def jacobi_update(self, int i, int j, int k):
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        return (1/6.) * (l1 + l2 + l3 + l4)

    def animate(self):
        anim = FuncAnimation(self.fig, self.update)
        plt.show()

    def run_animation(self):
        """
        Gives the ability to click on the animation canvas to play and pause.
        """
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
        anim = FuncAnimation(self.fig, self.update, interval=500)
        plt.show()
