import numpy as np
# cimport numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

# cdef class Grid:
class Grid:
    # cdef int N, M, steps_per_sweep
    # cdef float J, Kb, T
    # # cdef object anim , fig
    # cpdef object fig
    # cdef np.int64_t[:, :] grid
    # cpdef object __cpinit__(self, int N, int M, float T, float J=1, float Kb=1, all_up=True, anim=True):
    def __init__(self, int N, int M, float T, float J=1, float Kb=1, all_up=True, anim=True):
        self.N = N
        self.M = M
        self.J = J
        self.Kb = Kb
        self.T = T
        self.steps_per_sweep = N * M
        self.anim = anim
        if anim:
            self.fig = plt.figure()
        if all_up:
            # cdef np.int64_t[:, :]
            self.grid = np.ones((N, M))
        else:
            # cdef np.int64_t[:, :]
            self.grid = np.random.choice([-1, 1], size=(N, M))

    def print_grid(self):
        print(self.grid)

    def imshow_grid(self):
        plt.imshow(self.grid, interpolation='None', cmap='Blues', vmin=-1, vmax=1)

    def update_sweep(self, int k):
        cdef int N, M, i, n, m
        # N, M = self.grid.shape
        N = self.N
        M = self.M
        for i in range(self.steps_per_sweep):
            n = np.random.randint(N)
            m = np.random.randint(M)
            self.glauber_dynamics(n, m)
        if self.anim:
            self.fig.clear()
            self.imshow_grid()

    def glauber_dynamics(self, int n, int m):
        cdef int total, N, M, i, j
        # cdef float dE
        total = 0
        # N, M = self.grid.shape
        N = self.N
        M = self.M
        total += self.grid[(n - 1) % N][m]
        total += self.grid[(n + 1) % N][m]
        total += self.grid[n][(m-1) % M]
        total += self.grid[n][(m+1) % M]
        cdef float dE = 2 * self.J * self.grid[n][m] * total  # Check energy signs
        if dE <= 0:
            self.grid[n][m] *= -1
        elif np.random.rand() <= self.P(dE):
            self.grid[n][m] *= -1

    def P(self, float dE):
        cdef float exp
        if self.T == 0:
            return 0
        else:
            exp = np.exp (-dE / (self.Kb * self.T))
            return exp

    def temperature_tests(self, float t_min=1, float t_max=3, int data_points=20, int sweeps=100, int tests=50, eng=True, mag=True, save=True):
        cdef double [:] temperature = np.linspace(t_min, t_max, data_points)
        cdef double [:,:] magnetisation = np.zeros((data_points, tests))
        cdef double [:, :] energy = np.zeros((data_points, tests))
        cdef int i, j
        for i in range(data_points):
            sys.stdout.write("Simulation progress: %.1f%%\r" % (100 * i / data_points))
            sys.stdout.flush()

            self.T = temperature[i]  # Set the temperature of the system.
            self.update_sweep(sweeps)
            for j in range(tests):
                self.update_sweep(10)
                if mag:
                    magnetisation[i][j] = self.sys_magnetisation()
                if eng:
                    energy[i][j] = self.sys_energy()
        if save:
            np.savetxt('temperature.txt', temperature)
            if mag:
                np.savetxt('magnetisation.txt', magnetisation)
            if eng:
                np.savetxt('energy.txt', energy)

    def sys_magnetisation(self):
        """
        Magnetisation per data point in the simulation
        :return: <|m|> = |m_i/n|
        """
        cdef int M, mag
        M = np.sum(self.grid)
        mag = np.abs(M)
        return mag

    def sys_energy(self):
        # TODO needs to be tested
        cdef int N, M, n, m
        cdef float energy
        # N, M = self.grid.shape
        N = self.N
        M = self.M
        energy = 0
        for n in range(N):
            for m in range(M):
                ne_sum = self.grid[(n + 1) % N][m] + self.grid[n][(m + 1) % M]
                energy += -self.J * self.grid[n][m] * ne_sum  # Check if -J * ... or J * ... same with above in glauber
        return energy

    def susceptibility(self, save=True):
        cdef double [:,:] data
        cdef double [:] temp
        data = np.genfromtxt('magnetisation.txt')
        temp = np.genfromtxt('temperature.txt')
        # cdef double [:] magnetisation, chi
        cdef int x, i
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        chi = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N * self.M * self.Kb * temp[i])
            chi[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt('susceptibility.txt', chi)

    def heat_cap(self, save=True):
        cdef double [:, :] energy
        cdef double [:] C, temp
        cdef int i
        cdef double norm_fact
        data = np.genfromtxt('energy.txt')
        temp = np.genfromtxt('temperature.txt')
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        C = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N**2 * self.Kb * temp[i]**2)
            C[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt('heat_cap.txt', C)

    def animate(self):
        anim = FuncAnimation(self.fig, self.update_sweep)
        plt.show()
