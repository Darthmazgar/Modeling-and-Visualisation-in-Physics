"""
This provides a tiny change to the update methods where the laplacian is
not u(n+1) = u (n) + 1/dt**2 * (up+down+left+right-4*self) but is
u(n+1) = 1/dt**2 * (up+down+left+right-4*self)
This still dosent fix the issues.
"""

"""
2019 Modeling and Visualisation in Physics Exam
B075090

Update method not working but have attempted the rest.

There are two attempts at the update method. One using a rolling laplacian grid
and the other using point updates.

To change these change the update call in FuncAnimaiton between update and update_
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import exp
import sys

class Autocat:
    def __init__(self, N, R, F, D1=0.2, D2=0.1, k=0.06, dx=1, dt=0.01):
        self.N = N
        self.R = R
        self.F = F
        self.D1 = D1
        self.D2 = D2
        self.k = k
        self.dx = dx
        self.dt = dt
        self.u_grid = self.init_grid(greater_than=1, less_than=0.5)
        self.v_grid = self.init_grid(greater_than=0.01, less_than=0.25)
        self.animation = False
        self.sweeps_per_update = 1

    def init_grid(self, greater_than, less_than):
        """
        Initialises a grid such that one region less than self.R has the value
        less_than and other regions have the value greater_than.
        """
        grid = np.zeros((self.N, self.N))
        centre_x = int(self.N/2)
        centre_y = int(self.N/2)
        for i in range(self.N):
            for j in range(self.N):
                noise = np.random.uniform(low=-0.01, high=0.01)
                r = np.sqrt((i - centre_y)**2 + (j - centre_x)**2)
                if r >= self.R:
                    grid[i][j] = greater_than + noise
                else:
                    grid[i][j] = less_than + noise
        # Quick check to make sure that this works.
        # plt.imshow(grid)
        # plt.colorbar()
        # plt.show()
        return grid.copy()


    def laplacian_(self, grid, i, j):
        """
        Calculate the discrete laplacian for an individual point on a 2D grid
        with periodic boundary conditions.
        """
        l1 = grid[(i+1+self.N) % self.N][j] + grid[(i-1+self.N) % self.N][j]
        l2 = grid[i][(j+1+self.N) % self.N] + grid[i][(j-1+self.N) % self.N]
        l3 = -4*grid[i][j]
        return (l1 + l2 + l3)/self.dx**2

    def lap2D(self, lat):
        """
        Returns the (discretised) laplacian at each point in a 2d lattice with periodic BCs
        :param dx:   The size of finite steps in space
        :param lat:  The lattice of values in space for gradient calculation
        :return lap: Numpy array with values of the laplacian at each lattice point.
        """
        lap = np.roll(lat, 1, 0) + np.roll(lat, -1, 0) + \
              np.roll(lat, 1, 1) + np.roll(lat, -1, 1) - \
              4. * lat
        lap = 1./self.dx**2. * lap
        # print(lap[50][50])
        return(lap)

    def update(self, k):
        """
        Updates U and V grids by finding the change in each of the grids using
        a rolling laplacian grid method which updates the full grid in one
        operation.

        This appears to currently not be working hence the discrete update method below.

        Tried multiplying the full step by dt and this appears better but still wrong.

        I will move on and get the rest done before coming back to this.

        Have attempted the rest but without this update working there is no real
        way of checking if it is right.
        """
        for z in range(self.sweeps_per_update):


            u_update = self.dt* (np.multiply(self.D1,(self.lap2D(self.u_grid))))\
                    - self.dt * np.multiply(self.u_grid, np.square(self.v_grid))\
                     + self.dt * self.F * (1 - self.u_grid)


            v_update = self.dt*(np.multiply(self.D2,(self.lap2D(self.v_grid))))\
             + self.dt * np.multiply(self.u_grid, np.square(self.v_grid))\
              - self.dt * (self.F + self.k)*self.v_grid


            self.u_grid = np.add(self.u_grid, u_update)
            self.v_grid = np.add(self.u_grid, v_update)

        if self.animation:
            self.fig.clear()
            plt.xlabel("F: %.3f, dt:%.3f" % (self.F, self.dt))
            plt.imshow(self.u_grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower')
            plt.colorbar()

    def update_(self, k):
        """
        Alternative discrete update which I am still not getting working correctly at the moment.
        """
        for z in range(self.sweeps_per_update):
            new_u_grid = self.u_grid.copy()
            new_v_grid = self.v_grid.copy()
            for i in range(self.N):
                for j in range(self.N):

                    deltaU = (self.D1*self.dt) * (self.laplacian_(self.u_grid, i, j))\
                            - self.dt * self.u_grid[i][j]*self.v_grid[i][j]**2 \
                            + self.dt * self.F*(1-self.u_grid[i][j])
                    new_u_grid[i][j] += deltaU
                    deltaV = (self.D2*self.dt) * (self.laplacian_(self.v_grid, i, j))\
                            + self.dt*self.u_grid[i][j]*self.v_grid[i][j]**2 \
                            - self.dt*(self.F+self.k)*self.v_grid[i][j]
                    new_v_grid += deltaV
            self.u_grid = new_u_grid.copy()
            self.v_grid = new_v_grid.copy()
        if self.animation:
            self.fig.clear()
            plt.imshow(self.u_grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower')
            plt.colorbar()

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
                anim.event_source.stop()
                anim_running = False
            else:
                anim.event_source.start()
                anim_running = True

        self.fig.canvas.mpl_connect('button_press_event', onClick)
        anim = FuncAnimation(self.fig, self.update, fargs=None, interval=5)
        plt.show()

    def run_c_tests(self, n_tests=7, steps_per_measurement=10, measurements_per_test=500):
        self.dx = 1
        F_arr = np.linspace(0.020, 0.055, n_tests)
        var_results = np.zeros(n_tests)
        for i in range(n_tests):
            self.F  = F_arr[i]
            # Set up new grids for each test.
            self.u_grid = self.init_grid(greater_than=1, less_than=0.5)
            self.v_grid = self.init_grid(greater_than=0.01, less_than=0.25)

            # Initialising sweeps to allow for the sys to settle.
            self.sweeps_per_update=500
            self.update(1)
            self.sweeps_per_update = steps_per_measurement

            var = np.zeros(measurements_per_test)
            for j in range(measurements_per_test):
                var[j] = np.average(np.square(self.u_grid)) - np.square(np.average(self.u_grid))
                self.update(1)
            var_results[i] = np.average(var)
        var_results = np.array(list(zip(F_arr, var_results)))
        np.savetxt("var_results_c.dat", var_results, header="Varience data for part c.")
        plt.plot(var_results[:,0], var_results[:,1])
        plt.title("Varience changing with F")
        plt.xlabel("F,      k=%.3f" % self.k)
        plt.ylabel("Varience")
        plt.show()

    def run_d_tests(self, n_tests=7, steps_per_measurement=10, measurements_per_test=500):
        self.dx = 2
        self.k = 0.05
        F_arr = np.linspace(0.005, 0.030, n_tests)
        var_results = np.zeros(n_tests)
        error = np.zeros(n_tests)
        for i in range(n_tests):
            self.F = F_arr[i]
            # Set up new grids for each test.
            self.u_grid = self.init_grid(greater_than=1, less_than=0.5)
            self.v_grid = self.init_grid(greater_than=0.01, less_than=0.25)

            # Initialising sweeps to allow for the sys to settle.
            self.sweeps_per_update=500
            self.update(1)
            self.sweeps_per_update = steps_per_measurement

            var = np.zeros(measurements_per_test)
            for j in range(measurements_per_test):
                var[j] = np.average(np.square(self.u_grid)) - np.square(np.average(self.u_grid))
                self.update(1)
            var_results[i] = np.average(var)
            error[i] = self.bootstrap_errors(var)
        var_results = np.array(list(zip(F_arr, var_results, error)))
        np.savetxt("var_results_d.dat", var_results, header="Varience data for part d.\n[x,y,error]\nError calculated using bootstrap method.")
        plt.plot(var_results[:,0], var_results[:,1])
        plt.title("Varience changing with F")
        plt.xlabel("F,      k=%.3f" % self.k)
        plt.ylabel("Varience")
        plt.show()

    def varience(self, arr):
        return np.average(np.square(arr)) - np.square(np.average(arr))

    def bootstrap_errors(self, arr, k=100):
        """
        Calculate the errors of an array arr using the bootstrap error method
        of random sampling of the values to obtain effectivly more data.
        This works for one data point at a time and the values obtained for it.
        """
        val = np.zeros(k)

        for i in range(k):
            rand_selection = np.random.choice(arr, len(arr))
            val[i] = self.varience(rand_selection)
        sigma = np.sqrt(np.average(np.square(val))
                 - np.square(np.average(val)))
        return sigma
