import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
from libcpp cimport bool
from scipy.ndimage import generate_binary_structure

class Grid:
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
            self.grid = np.ones((N, M))
        if not all_up:
            self.grid = np.random.choice([-1, 1], size=(N, M))

    def print_grid(self):
        print(self.grid)

    def imshow_grid(self):
        plt.imshow(self.grid, interpolation='None', cmap='Blues', vmin=-1, vmax=1)

    def update_sweep(self, int k):
        cdef int N, M, i, n, m
        N, M = self.grid.shape
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
        N, M = self.grid.shape
        k = generate_binary_structure(2, 1)
        for i in range(n-1, n+2):
            for j in range(m-1, m+2):
                # if k[i - n + 1][j - m + 1]:
                #     total += self.grid[i % N][j % M]
                if i != n or j != m:
                    # continue
                    total += self.grid[i % N][j % M]
                # if i == n or j == n:
                #     # print(i, n)
                #     # print(j, m)
                #     total += self.grid[i % N][j % M]
        # print("total: %d" % total)
        cdef float dE = 2 * self.J * self.grid[n][m] * total  # Check energy signs
        if dE <= 0:
            self.grid[n][m] *= -1
        elif np.random.rand() <= self.P(dE):
            # print(np.exp(-dE / (self.Kb * self.T)))
            self.grid[n][m] *= -1

    def P(self, float dE):
        cdef float exp
        if self.T == 0:
            # print("T=0")
            return 0
        else:
            exp = np.exp (-dE / (self.Kb * self.T))
            # print(exp)
            return exp

    def temperature_tests(self, float t_min=2, float t_max=10, int data_points=10, int sweeps=50, int tests=10, eng=True, mag=True, save=True):
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
        N, M = self.grid.shape
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

# def main():
#     grid = Grid(50, 50, 0, anim=True, all_up=False)
#     # grid.print_grid()
#     # grid.imshow_grid()
#     # plt.show()
#     for i in range(100):
#         grid.update_sweep(1)
#     grid.imshow_grid()
#     # plt.show()
#     # grid.animate()
#     # grid.temperature_tests()
#     # grid.susceptibility()
#     # grid.heat_cap()

# main()
