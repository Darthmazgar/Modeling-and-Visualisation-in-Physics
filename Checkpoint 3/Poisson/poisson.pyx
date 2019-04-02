import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys

# class Vector:
#     def __init__(self, x, y, z):
#         self.x = x
#         self.y = y
#         self.z = z

class Poisson:
    def __init__(self, int N, float accuracy, float dx=1, float epsilon=1):
        self.N = N
        self.acc = accuracy
        self.dx = dx
        self.epsilon = epsilon
        self.rho_grid = np.zeros((N, N, N))
        self.next_rho_grid = np.zeros((N, N, N))
        self.phi_grid = np.zeros((N, N, N))
        self.next_phi_grid = np.zeros((N, N, N))

    def zero_boundaries(self):
        self.phi_grid[:, :, 0] = 0
        self.phi_grid[:, :, -1] = 0
        self.phi_grid[:, 0, :] = 0
        self.phi_grid[:, -1, :] = 0
        self.phi_grid[0, :, :] = 0
        self.phi_grid[-1, :, :] = 0

    def add_point_charge(self, int i, int j, int k):
        if i >= self.N or j >= self.N or k >= self.N:
            raise ValueError("Point selected (%d, %d, %d) out of range" % (i, j, k)
                        + " (%d, %d, %d)" % (self.N, self.N, self.N))
        self.rho_grid[i][j][k] = 1

    def add_line_charge(self, int i, int j):
        if i >= self.N or j >= self.N:
            raise ValueError("Point selected (%d, %d) out of range" % (i, j)
                        + " (%d, %d)" % (self.N, self.N))
        for z in range(1, self.N-1):  # Maybe from 0 -> self.N
            self.rho_grid[z][i][j] = 1

    def update(self, int k):
        for z in range(1000):  # sweeps
            for i in range(1, self.N-1):
                for j in range(1, self.N-1):  # 1 -> N-1 to preserve zero at boundary
                    for k in range(1, self.N-1):
                        # self.next_phi_grid[i][j][k] = self.jacobi_update(i, j, k)
                        self.next_phi_grid[i][j][k] = self.gauss_seidel_update(i, j, k)

            self.phi_grid = self.next_phi_grid.copy()

        # self.fig.clear()
        # plt.imshow(self.phi_grid[12], interpolation='nearest',
        #                cmap='coolwarm', origin='lower')
        # plt.colorbar()

    def jacobi_update(self, int i, int j, int k):
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        # val = (1/6.) * (l1 + l2 + l3 + l4)
        # if val != 0:
        #     print(val)
        return (1/6.) * (l1 + l2 + l3 + l4)

    def gauss_seidel_update(self, int i, int j, int k):
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.next_phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.next_phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.next_phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        return (1/6.) * (l1 + l2 + l3 + l4)

    def rho_update(self, int i, int j, int k):
        cdef double grad_sq_phi
        grad_sq_phi = (self.phi_grid[(i+1 + self.N) % self.N][j][k]
            + self.phi_grid[(i-1 + self.N) % self.N][j][k]
            + self.phi_grid[i][(j+1 + self.N) % self.N][k]
            + self.phi_grid[i][(j-1 + self.N) % self.N][k]
            + self.phi_grid[i][j][(k+1 + self.N) % self.N]
            + self.phi_grid[i][j][(k-1 + self.N) % self.N]
            - 6 * self.phi_grid[i][j][k])
        return (- grad_sq_phi * self.epsilon)

    # def e_field(self, int i, int j, int k):  # Made obselete by using np.gradient(self.phi_grid)
    #     dx = self.phi_grid[(i+1+self.N) % self.N][j][k] - self.phi_grid[(i-1+self.N) % self.N][j][k]
    #     dy = self.phi_grid[i][(j+1+self.N) % self.N][k] - self.phi_grid[i][(j-1+self.N) % self.N][k]
    #     dz = self.phi_grid[i][j][(k+1+self.N) % self.N] - self.phi_grid[i][j][(k-1+self.N) % self.N]
    #     return np.array((dx, dy, dz))

    def plot_E_field(self):
        E_field = np.gradient(self.phi_grid)
        # print(np.shape(E_field))
        ex = E_field[2][int(self.N/2)][:][:]
        ey = E_field[1][int(self.N/2)][:][:]
        plt.quiver(ex, ey, pivot='tip')
        plt.show()

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
