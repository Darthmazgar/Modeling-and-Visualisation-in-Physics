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
    def __init__(self, int N, float accuracy, method, float dx=1, float epsilon=1):
        self.N = N
        self.acc = accuracy
        self.dx = dx
        self.epsilon = epsilon
        self.sweeps = 1
        self.method = method
        self.rho_grid = np.zeros((N, N, N))
        self.next_rho_grid = np.zeros((N, N, N))
        self.phi_grid = np.zeros((N, N, N))
        self.next_phi_grid = np.zeros((N, N, N))

        self.animation = False

    def zero_boundaries(self):
        self.phi_grid[:, :, 0] = 0
        self.phi_grid[:, :, -1] = 0
        self.phi_grid[:, 0, :] = 0
        self.phi_grid[:, -1, :] = 0
        self.phi_grid[0, :, :] = 0
        self.phi_grid[-1, :, :] = 0

    def add_point_charge(self, int i, int j, int k, float val):
        if i >= self.N or j >= self.N or k >= self.N:
            raise ValueError("Point selected (%d, %d, %d) out of range" % (i, j, k)
                        + " (%d, %d, %d)" % (self.N, self.N, self.N))
        self.rho_grid[i][j][k] = val

    def add_line_charge(self, int i, int j, float val):
        if i >= self.N or j >= self.N:
            raise ValueError("Point selected (%d, %d) out of range" % (i, j)
                        + " (%d, %d)" % (self.N, self.N))
        for z in range(1, self.N-1):  # Maybe from 0 -> self.N
            self.rho_grid[z][i][j] = val
        # self.rho_grid[i][i][:] = val

    def update(self, int k):
        cdef int z, i, j
        for z in range(self.sweeps):  # sweeps
            for i in range(1, self.N-1):
                for j in range(1, self.N-1):  # 1 -> N-1 to preserve zero at boundary
                    for k in range(1, self.N-1):
                        if self.method == 'jacobi':
                            self.next_phi_grid[i][j][k] = self.jacobi_update(i, j, k)
                        else:
                            self.next_phi_grid[i][j][k] = self.gauss_seidel_update(i, j, k)

            self.phi_grid = self.next_phi_grid.copy()

        if self.animation:
            self.fig.clear()
            plt.imshow(self.phi_grid[int(self.N/2)], interpolation='nearest',
                           cmap='coolwarm', origin='lower')
            plt.colorbar()

    def jacobi_update(self, int i, int j, int k):
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        return (1/6.) * (l1 + l2 + l3 + l4)

    def gauss_seidel_update(self, int i, int j, int k):
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.next_phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.next_phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.next_phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        return (1/6.) * (l1 + l2 + l3 + l4)

    def plot_E_field(self, save=False):
        E_field = np.gradient(self.phi_grid)
        ex = E_field[2][int(self.N/2)][:][:]
        ey = E_field[1][int(self.N/2)][:][:]
        if save:
            np.savetxt('ex.txt', ex, header=('The x component of the electric field.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps)))
            np.savetxt('ey.txt', ex, header=('The y component of the electric field.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps)))
        plt.quiver(ex, ey, pivot='tip')
        plt.title('Electric field')
        plt.xlabel("Grid size: {}$^3$, Init sweeps: {}".format(self.N, self.sweeps))
        plt.show()

    def plot_B_field(self, save=False):
        if np.count_nonzero(self.rho_grid) == 1:
            raise ValueError("No magnetic monopoles exist. Try adding more charges.")
        bx, by = np.gradient(self.phi_grid[int(self.N/2)][:][:], edge_order=0)
        if save:
            np.savetxt('bx.txt', bx, header=('The x component of the magnetic field.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps)))
            np.savetxt('by.txt', by, header=('The y component of the magnetic field.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps)))

        plt.quiver(bx, -by, pivot='tip')
        plt.title('Magnetic field')
        plt.xlabel("Grid size: {}$^3$, Init sweeps: {}".format(self.N, self.sweeps))
        plt.show()

    def animate(self):
        anim = FuncAnimation(self.fig, self.update)
        plt.show()

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
