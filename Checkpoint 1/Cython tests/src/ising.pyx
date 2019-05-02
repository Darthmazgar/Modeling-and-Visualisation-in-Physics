import numpy as np
# cimport numpy as np
import cython
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
from libc.stdlib cimport abs, rand
from libc.math cimport exp

class IsingGrid:
    def __init__(self, N, M, T, ds=ds, sv_ext='', anim=True, all_up=False):
        self.N = N
        self.M = M
        self.T = T
        self.ds = ds
        self.sv_ext = sv_ext
        self.anim = anim
        self.all_up = all_up


        
        plt.imshow(self.grid, interpolation='sinc',
                   cmap='Blues', vmin=-1, vmax=1)

    def update_sweep(self, int k):
        """
        Runs through one sweep of changes in the grid.
        :param: (int) k, variable used to count frames by FuncAnimation.
        """
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

    def nn_check(self, int n, int m, int x, int y):
        """
        Check if two points (x,y) and (n,m) are nearest neighbours.
        :return: (Boolean) True if nearest neighbours; False otherwise.
        """
        cdef int N, M
        N = self.N
        M = self.M
        if n == x and (m == y or m == (y+1) % N or m == (y-1) % M):
            return True
        elif m == y and (n == (x+1) % N or n == (x-1) % M):
            return True
        else:
            return False

    def kawasaki_dynamics(self):
        """
        Updates the grid according to Kawasaki synamics. Two points are
        considred randomly. If it is energeticaly favourable for those points to
        be switched then they are. If not then there is still a probability that
        they are switched which varies with the temperature of the
        system (self.T).
        :return: 1 if swaped; 0 if left the same.
        """
        cdef int total, total2, N, M, x, y
        cdef float dE
        N = self.N
        M = self.M
        n = rand() % N
        m = rand() % M
        x = rand() % N
        y = rand() % M
        if self.grid[n][m] == self.grid[x][y]:
            return 0
            # If the points are the same then changing would do nothing.
        total = self.sum_spin(n, m, N, M)
        total2 = self.sum_spin(x, y, N, M)
        if self.nn_check(n, m, x, y):
            # If the points are nearest neighbours then the 'total' of the
            # points around is increased by 4 to account for the double
            # counting.
            dE = 2 * self.J * (self.grid[n][m]*total + self.grid[x][y]*total2)\
                + 4* self.J
        else:
            # Calculate the energy change
            dE = 2 * self.J * (self.grid[n][m]*total + self.grid[x][y]*total2)
        if dE <= 0:  # Swap spin of both if favourable.
            self.grid[n][m] *= -1
            self.grid[x][y] *= -1
        elif np.random.rand() <= self.P(dE):
            # Swap spin of both if random is less than probability from
            # self.P().
            self.grid[n][m] *= -1
            self.grid[x][y] *= -1
        return 1


    def glauber_dynamics(self):
        """
        Flipes the state of a randomly selected point if it is either
        energeticaly favourable to do so or based on a probability which depends
        on the temperature of the system.
        :return: 1 if swaped; 0 if left the same.
        """
        cdef int total, N, M
        N = self.N
        M = self.M
        n = rand() % N
        m = rand() % M
        total = self.sum_spin(n, m, N, M)  # Sum values of points around.
        cdef float dE = 2 * self.J * self.grid[n][m] * total
        if dE <= 0:  # Flip if energeticaly favourable.
            self.grid[n][m] *= -1
            return 1
        elif np.random.rand() <= self.P(dE):
            # Flip if rand is less than the probability given by self.P().
            self.grid[n][m] *= -1
            return 1
        else:
            return 0

    def sum_spin(self, int n, int m, int N, int M):
        """
        Sums the values of the points in the grid around a selected point (n, m)
        in a grid size of N*M.
        :return: The sum of the values.
        """
        cdef int total = 0
        total += self.grid[(n - 1) % N][m]
        total += self.grid[(n + 1) % N][m]
        total += self.grid[n][(m - 1) % M]
        total += self.grid[n][(m + 1) % M]
        return total

    @cython.cdivision(True)
    def P(self, float dE):
        cdef float expo
        if self.T == 0:
            return 0
        else:
            expo = exp(-dE / (self.Kb * self.T))
            return expo

    def temperature_tests(self, float t_min=1, float t_max=2.9, int data_points=20,
                        int sweeps=100, int tests=10000, int sweeps_per_test=10,
                        eng=True, mag=True, save=True):

        cdef double [:] temperature = np.linspace(t_min, t_max, data_points)
        cdef double [:,:] magnetisation = np.zeros((data_points, tests))
        cdef double [:, :] energy = np.zeros((data_points, tests))
        cdef int i, j
        for i in range(data_points):
            sys.stdout.write("Simulation progress: %.1f%%\r"
                            % (100 * i / data_points))
            sys.stdout.flush()  # Prints progress of simulation.

            self.T = temperature[i]  # Set the temperature of the system.
            self.update_sweep(sweeps * 2)  # Run simulation for a given number
                                           # of sweeps. x2 longer init time.
            for j in range(tests):
                self.update_sweep(sweeps_per_test)
                if mag:
                    # Calculate current system magnetisation.
                    magnetisation[i][j] = self.sys_magnetisation()
                if eng:
                    # Calculate current system energy.
                    energy[i][j] = self.sys_energy()
        if save:
            np.savetxt(('temperature'+ self.sv_ext +'.txt'), temperature)
            if mag:
                np.savetxt(('magnetisation'+ self.sv_ext +'.txt'),magnetisation)
            if eng:
                np.savetxt(('energy'+ self.sv_ext +'.txt'), energy)

    def sys_magnetisation(self):
        """
        Magnetisation per data point in the simulation
        :return: <|m|> = |m_i/n|
        """
        cdef int M, mag
        M = np.sum(self.grid)
        mag = abs(M)
        return mag

    def sys_energy(self):
        """
        Calculate the system energy at the current state.
        :return: (float) System energy.
        """
        cdef int N, M, n, m
        cdef float energy
        N = self.N
        M = self.M
        energy = 0
        for n in range(N):
            for m in range(M):
                ne_sum = self.grid[(n + 1) % N][m] + self.grid[n][(m + 1) % M]
                energy += -self.J * self.grid[n][m] * ne_sum
        return energy

    def susceptibility(self, save=True):
        """
        Reads in the magnetisation and temperature from saved text files and
        calculates the susceptibility at each temperature.
        :param: (Boolean) save: saves the susceptibility measurments to a text
        file if True.
        """
        cdef double [:,:] data
        cdef double [:] temp
        data = np.genfromtxt(('magnetisation'+ self.sv_ext +'.txt'))
        temp = np.genfromtxt(('temperature'+ self.sv_ext +'.txt'))

        cdef int x, i
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        chi = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N * self.M * self.Kb * temp[i])
            chi[i] = norm_fact * (np.average(np.square(data[i]))
                   - np.square(np.average(data[i])))
        if save:
            np.savetxt(('susceptibility'+ self.sv_ext +'.txt'), chi)
        return chi

    def heat_cap(self, save=True):
        """
        Reads in the energy and temperature from saved text files and
        calculates the heat capacity at each temperature.
        :param: (Boolean) save: saves the heat capacity measurments to a text
        file if True.
        """
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
            C[i] = norm_fact * (np.average(np.square(data[i]))
                 - np.square(np.average(data[i])))
        if save:
            np.savetxt(('heat_cap'+ self.sv_ext +'.txt'), C)
        return C

    def bootstarap_errors(self, int k=100, save=True):
        """
        Reads in the energy and temperature from saved text files and
        calculates the ascosiated error at each temperature according to the
        'bootstrap' method for calculating errors. N points are chosen at random
        from each energy reading. This is done 'k' times. These results are used
        to find the varience at each step to give an error.
        :param: (Boolean) save: saves the ascosiated error measurments to a text
        file if True.
        """
        cdef int i, dlen
        cdef double avg, norm_fact
        cdef double [:, :] data
        cdef double [:] heat_cap, temp ,sel_data, h_caps, sigma
        data = np.genfromtxt(('energy'+ self.sv_ext +'.txt'))
        temp = np.genfromtxt(('temperature'+ self.sv_ext +'.txt'))
        dlen = len(data)
        row_len = len(data[0])
        heat_cap = np.zeros(k)
        sigma = np.zeros(dlen)
        for i in range(dlen):
            norm_fact = 1 / (self.N**2 * self.Kb * temp[i]**2)
            for j in range(k):
                sel_data = np.random.choice(data[i], row_len)
                # Choose N rand points.
                heat_cap[j] = norm_fact * (np.average(np.square(sel_data))
                              - np.square(np.average(sel_data)))
            sigma[i] = np.sqrt(np.average(np.square(heat_cap))
                     - np.square(np.average(heat_cap)))
        if save:
            np.savetxt(('sigma_bs' + self.sv_ext + '.txt'), sigma)
        return sigma

    def animate(self):
        anim = FuncAnimation(self.fig, self.update_sweep)
        plt.show()
