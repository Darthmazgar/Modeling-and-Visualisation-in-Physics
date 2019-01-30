import numpy as np
# cimport numpy as np
import cython
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
import sys
from libc.stdlib cimport abs, rand
from libc.math cimport exp

# cdef class Grid:
class Grid:
    # cdef int N, M, steps_per_sweep
    # cdef float J, Kb, T
    # # cdef object anim , fig
    # cpdef object fig
    # cdef np.int64_t[:, :] grid
    # cpdef object __cpinit__(self, int N, int M, float T, float J=1, float Kb=1, all_up=True, anim=True):
    def __init__(self, int N, int M, float T, sv_ext, ds=0, float J=1, float Kb=1, all_up=True, anim=True):
        self.N = N
        self.M = M
        self.J = J
        self.Kb = Kb
        self.T = T
        self.ds = ds
        self.sv_ext = sv_ext
        self.steps_per_sweep = N * M
        self.anim = anim
        if anim:
            self.fig = plt.figure()
        if all_up:
            self.grid = np.ones((N, M))
        else:
            self.grid = np.random.choice([-1, 1], size=(N, M))

    def print_grid(self):
        print(self.grid)

    def imshow_grid(self):
        # axcolor = 'lightgoldenrodyellow'
        # axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
        plt.imshow(self.grid, interpolation='None', cmap='Blues', vmin=-1, vmax=1)
        # self.T = Slider(axfreq, 'Temperature', 0, 4, valinit=0.8)

    def update_sweep(self, int k):
        cdef int N, M, i, n, m
        N = self.N
        M = self.M
        for i in range(self.steps_per_sweep):
            if self.ds == 0:
                self.glauber_dynamics()
            elif self.ds == 1:
                self.kawasaki_dynamics()
        if self.anim:
            self.fig.clear()
            self.imshow_grid()

    def nn_check(self, n, m, x, y):
        if n == x and (m == y or m == y+1 or m == y-1):
            return True
        elif m == y and (n == x+1 or n == x-1):
            return True
        else:
            return False

    def kawasaki_dynamics(self):
        cdef int total, total2, N, M, x, y
        cdef float dE
        # n = np.random.randint(N)
        # m = np.random.randint(M)
        N = self.N
        M = self.M
        n = rand() % N
        m = rand() % M
        x = rand() % N
        y = rand() % M
        if self.grid[n][m] == self.grid[x][y]:
            return 0
        total = self.sum_spin(n, m, N, M)
        total2 = self.sum_spin(x, y, N, M)
        # cdef float ei = -self.J * self.grid[n][m] * total - self.J * self.grid[x][y] * total2
        # cdef float ef = self.J * self.grid[n][m] * total + self.J * self.grid[x][y] * total2
        if self.nn_check(n, m, x, y):
            # TODO Need to work out what goes on with nn
            dE = 0 # ###############
        else:
            dE = 2 * self.J * (self.grid[n][m]*total + self.grid[x][y]*total2)
        # print(dE)
        if dE <= 0:
            self.grid[n][m] *= -1
            self.grid[x][y] *= -1
        elif np.random.rand() <= self.P(dE):
            self.grid[n][m] *= -1
            self.grid[x][y] *= -1
        return 1


    def glauber_dynamics(self):
        cdef int total, N, M
        # n = np.random.randint(N)
        # m = np.random.randint(M)
        N = self.N
        M = self.M
        n = rand() % N
        m = rand() % M
        total = self.sum_spin(n, m, N, M)
        cdef float dE = 2 * self.J * self.grid[n][m] * total  # Check energy signs
        if dE <= 0:
            self.grid[n][m] *= -1
        elif np.random.rand() <= self.P(dE):
            self.grid[n][m] *= -1

    def sum_spin(self, int n, int m, int N, int M):
        cdef int total = 0
        total += self.grid[(n - 1) % N][m]
        total += self.grid[(n + 1) % N][m]
        total += self.grid[n][(m-1) % M]
        total += self.grid[n][(m+1) % M]
        return total

    @cython.cdivision(True)
    def P(self, float dE):
        cdef float expo
        if self.T == 0:
            return 0
        else:
            expo = exp(-dE / (self.Kb * self.T))
            return expo

    def temperature_tests(self, float t_min=1, float t_max=3, int data_points=20,
                        int sweeps=100, int tests=1000, int sweeps_per_test=10,
                        eng=True, mag=True, save=True):
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
            np.savetxt(('temperature'+ self.sv_ext +'.txt'), temperature)
            if mag:
                np.savetxt(('magnetisation'+ self.sv_ext +'.txt'), magnetisation)
            if eng:
                np.savetxt(('energy'+ self.sv_ext +'.txt'), energy)

    def sys_magnetisation(self):
        """
        Magnetisation per data point in the simulation
        :return: <|m|> = |m_i/n|
        """
        cdef int M, mag
        M = np.sum(self.grid)
        # mag = np.abs(M)
        mag = abs(M)
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
        data = np.genfromtxt(('magnetisation'+ self.sv_ext +'.txt'))
        temp = np.genfromtxt(('temperature'+ self.sv_ext +'.txt'))
        # cdef double [:] magnetisation, chi
        cdef int x, i
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        chi = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N * self.M * self.Kb * temp[i])
            chi[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt(('susceptibility'+ self.sv_ext +'.txt'), chi)

    def heat_cap(self, save=True):
        cdef double [:, :] energy
        cdef double [:] C, temp
        cdef int i
        cdef double norm_fact
        data = np.genfromtxt(('energy'+ self.sv_ext +'.txt'))
        temp = np.genfromtxt(('temperature'+ self.sv_ext +'.txt'))
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        C = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N**2 * self.Kb * temp[i]**2)
            C[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt(('heat_cap'+ self.sv_ext +'.txt'), C)

    def animate(self):
        anim = FuncAnimation(self.fig, self.update_sweep)
        plt.show()
