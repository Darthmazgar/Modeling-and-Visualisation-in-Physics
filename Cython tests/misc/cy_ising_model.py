import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class Grid:
    def __init__(self, N, M, T, J=1, Kb=1, all_up=True, anim=True):
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

    def update_sweep(self, k):
        N, M = self.grid.shape
        for i in range(self.steps_per_sweep):
            n = np.random.randint(N)
            m = np.random.randint(M)
            self.glauber_dynamics(n, m)
        if self.anim:
            self.fig.clear()
            self.imshow_grid()

    def glauber_dynamics(self, n, m):
        # TODO i think there is something wrong here that makes the need for a much higher T
        total = 0
        N, M = self.grid.shape
        for i in range(n-1, n+2):
            for j in range(m-1, m+2):
                if i == n and j == m:
                    continue
                total += self.grid[i % N][j % M]
        dE = 2 * self.J * self.grid[n][m] * total  # Check energy signs
        if dE <= 0:
            self.grid[n][m] *= -1
        elif np.random.rand() <= self.P(dE):
            # print(np.exp(-dE / (self.Kb * self.T)))
            self.grid[n][m] *= -1

    def P(self, dE):
        if self.T == 0:
            return 0
        return np.exp(- (dE / (self.Kb * self.T)))

    def temperature_tests(self, t_min=1, t_max=3, data_points=2, sweeps=10, tests=1, eng=True, mag=True, save=True):
        temperature = np.linspace(t_min, t_max, data_points)
        magnetisation = np.zeros((data_points, tests))
        energy = np.zeros((data_points, tests))

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
        M = np.sum(self.grid)
        mag = np.abs(M)
        return mag

    def sys_energy(self):
        # TODO needs to be tested
        N, M = self.grid.shape
        energy = 0
        for n in range(N):
            for m in range(M):
                ne_sum = self.grid[(n + 1) % N][m] + self.grid[n][(m + 1) % M]
                energy += -self.J * self.grid[n][m] * ne_sum
        print(energy)
        return energy

    def susceptibility(self, save=True):
        data = np.genfromtxt('magnetisation.txt')
        temp = np.genfromtxt('temperature.txt')
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        chi = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N * self.M * self.Kb * temp[i])
            chi[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt('susceptibility.txt', chi)

    def heat_cap(self, save=True):
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

def main():
    grid = Grid(50, 50, 3, anim=True, all_up=True)
    # grid.print_grid()
    # grid.imshow_grid()
    # plt.show()
    for i in range(100):
        grid.update_sweep(1)
    grid.imshow_grid()
    # plt.show()
    # grid.animate()
    # grid.temperature_tests(save=False)
    # grid.susceptibility()
    # grid.heat_cap()

# main()
