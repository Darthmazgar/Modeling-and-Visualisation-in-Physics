import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys


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
        self.lim_reached = False

    def zero_boundaries(self):
        """ Sets the boundaries of the volume to be zero. """
        self.phi_grid[:, :, 0] = 0
        self.phi_grid[:, :, -1] = 0
        self.phi_grid[:, 0, :] = 0
        self.phi_grid[:, -1, :] = 0
        self.phi_grid[0, :, :] = 0
        self.phi_grid[-1, :, :] = 0

    def add_point_charge(self, int i, int j, int k, float val):
        """ Adds a point charge at point (i, j, k) where i, j and k are integers. """
        if i >= self.N or j >= self.N or k >= self.N:
            raise ValueError("Point selected (%d, %d, %d) out of range" % (i, j, k)
                        + " (%d, %d, %d)" % (self.N, self.N, self.N))
        self.rho_grid[i][j][k] = val

    def add_line_charge(self, int i, int j, float val):
        """ Adds a line of charge in aligned with the z-axis at position (i, j). """
        if i >= self.N or j >= self.N:
            raise ValueError("Point selected (%d, %d) out of range" % (i, j)
                        + " (%d, %d)" % (self.N, self.N))
        for z in range(1, self.N-1):  # Maybe from 0 -> self.N
            self.rho_grid[z][i][j] = val

    def add_plane(self, int y, float val):
        """ Adds a plane of charge in the x,y plane. """
        for z in range(self.N):
            for x in range(self.N):
                self.rho_grid[x][y][z] = val

    def add_gausian(self, float size, float val):
        """
        Adds a Gaussian distribution of charge in the centre of the space with
        size being a flot between 0 and 1 definin the sizr of the distribution
        for 1 sigma.
        """
        x, y, z = np.meshgrid(np.linspace(-1, 1, self.N), np.linspace(-1, 1, self.N), np.linspace(-1, 1, self.N))
        d = np.sqrt(x**2 + y**2 + z**2)
        sigma, mu = size, 0.1
        # print(np.shape(np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )))
        # print(np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) ))
        self.rho_grid = np.exp(-((d - mu)**2 / ( 2.0*sigma**2))) * val

    def update(self, int k):
        """
        Updates the phi grid for the full volume.
        :param: <int> k: Needed as frame number for animation. Set to 1 else.
        """
        cdef int z, i, j
        for z in range(self.sweeps):  # sweeps
            for i in range(1, self.N-1):
                for j in range(1, self.N-1):  # 1 -> N-1 to preserve zero at boundary
                    for k in range(1, self.N-1):
            # for i in range(self.N):
            #     for j in range(self.N):
            #         for k in range(self.N):
                        if self.method == 'jacobi':
                            self.next_phi_grid[i][j][k] = self.jacobi_update(i, j, k)
                        else:
                            self.next_phi_grid[i][j][k] = self.gauss_seidel_update(i, j, k)
            max_diff = np.max(np.sqrt(np.square(np.subtract(self.phi_grid, self.next_phi_grid))))
            if max_diff <= self.acc and not self.lim_reached:
                if not self.animation:
                    print("Accuracy level of {:0.6f} reached after {} sweeps.".format(self.acc, (z+1)))
                else:
                    print("Accuracy level of {:0.6f}.".format(self.acc))
                self.lim_reached = True
                break
            self.phi_grid = self.next_phi_grid.copy()

        if self.animation:  # If animating then plot slice.
            self.fig.clear()
            plt.title("Potential Strength")
            plt.imshow(self.phi_grid[int(self.N/2)], interpolation='nearest',
                           cmap='coolwarm', origin='lower')
            plt.colorbar()


    def jacobi_update(self, int i, int j, int k):
        """ Updates the point (i, j, k) following the Jacobi algorithm. """
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        return (1/6.) * (l1 + l2 + l3 + l4)

    def gauss_seidel_update(self, int i, int j, int k):
        """ Updates the point (i, j, k) following the Gauss-Seidel algorithm. """
        cdef double l1, l2, l3, l4
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j][k] + self.next_phi_grid[(i - 1 + self.N) % self.N][j][k]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N][k] + self.next_phi_grid[i][(j - 1 + self.N) % self.N][k]
        l3 = self.phi_grid[i][j][(k + 1 + self.N) % self.N] + self.next_phi_grid[i][j][(k - 1 + self.N) % self.N]
        l4 = self.dx**2 * self.rho_grid[i][j][k]
        return (1/6.) * (l1 + l2 + l3 + l4)

    def contour_test(self):
        """
        Produces a contour plot for the potential strength in space, along with
        a slice plot through the mid plane which can be used to check the fall
        of rate of the field strength.
        """
        self.sweeps = 100
        self.update(1)
        plt.imshow(self.phi_grid[int(self.N/2)], interpolation='nearest',
                       cmap='coolwarm', origin='lower')
        plt.colorbar()
        # print(np.shape(self.phi_grid[int(self.N/2)]))
        np.savetxt('point_contour.txt', self.phi_grid[int(self.N/2)], header='Point charge contour data.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps))
        plt.show()

        cut = self.phi_grid[int(self.N/2)][int(self.N/2)]
        xs = np.linspace(-int(self.N/2), int(self.N/2), self.N)
        np.savetxt('point_cut.txt', np.array(list(zip(xs,cut))), header='Point charge midplane cut data.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps))
        plt.plot(xs, cut)
        plt.title("Cut through point charge 1D")
        plt.xlabel("Displacement")
        plt.ylabel("Charge")
        plt.show()


    def plot_E_field(self, save=False):
        """ Plots the current electric field of the system as a slice through
        the x,y plane. """
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
        """ Plots the current magnetic field of the system as a slice through
        the x,y plane. If there is only one charge present reises a ValueError
        as no magnetic monopoles exist. """
        if np.count_nonzero(self.rho_grid) == 1:
            raise ValueError("No magnetic monopoles exist. Try adding more charges.")
        bx, by = np.gradient(self.phi_grid[int(self.N/2)][:][:], edge_order=0)
        if save:
            np.savetxt('bx.txt', bx, header=('The x component of the magnetic field.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps)))
            np.savetxt('by.txt', by, header=('The y component of the magnetic field.\nFor Grid size: {} and Init sweeps: {}'.format(self.N, self.sweeps)))
        if self.animation:
            plt.imshow(self.phi_grid[int(self.N/2)], interpolation='nearest',
                           cmap='coolwarm', origin='lower')
            plt.colorbar()
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
