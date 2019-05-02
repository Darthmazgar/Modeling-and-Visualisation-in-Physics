import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import exp
import sys

class SuperfluidTransiton: ################
    def __init__(self, N, T, c, J=-1, kb=1):
        self.N = N
        self.T = T
        self.c = c
        self.J = J
        self.kb = kb
        n_spins = 50
        spins = [x for x in range(-n_spins, n_spins)]
        self.grid = np.random.choice(spins, size=(N,N) )# , p=[c*0.5,c*0.5,1-c])
        self.animation = False
        self.sweeps_per_update = 1

    def dE(self, i, j, home=None, swap=None):
        """ Calculates the vhange in Energy of a given configuration by summing all
        nearest neighbours multiplied by the home cell value. """
        sum = 0
        if home == None:
            home = self.grid[i][j]
        if swap == None:
            swap = -1*home
        sum += self.grid[(i+1+self.N)%self.N][j]# *home
        sum += self.grid[(i-1+self.N)%self.N][j]# *home
        sum += self.grid[i][(j+1+self.N)%self.N]# *home
        sum += self.grid[i][(j-1+self.N)%self.N]# *home
        return -self.J * sum * (home - swap)

    def P(self, dE):
        """ Calculates the probability of being in a given configuration. """
        if self.T == 0:
            return 0
        else:
            return exp(-dE/(self.kb*self.T))

    def selection_update(self):
        He4 = False
        # ensure a He4 point is selected.
        while not He4:
            i, j = np.random.randint(self.N, size=2)
            val = self.grid[i][j]
            if val != 0:
                He4 = True
        dE = self.dE(i, j, val)
        if dE <= 0 or self.P(dE) > np.random.uniform():
            self.grid[i][j] *= -1
            return True
        else:
            return False

    def pair_selection_update(self):
        i, j, x, y = np.random.randint(self.N, size=4)
        val1 = self.grid[i][j]
        val2 = self.grid[x][y]
        dE1 = self.dE(i, j, val1, val2)
        dE2 = self.dE(x, y, val2, val1)
        dE = dE1 + dE2
        if dE <= 0 or self.P(dE) > np.random.uniform():
            self.grid[i][j] = val2
            self.grid[x][y] = val1
            return True
        else:
            return False

    def update(self, k):
        for z in range(self.sweeps_per_update):
            update_rule = np.random.choice([0, 1], size=self.N**2)
            # update_rule = [1 for x in range(self.N**2)]  # Set all updates to the same rule for testing.
            for i in range(self.N**2):
                if update_rule[i] == 0:
                    self.selection_update()
                elif update_rule[i] == 1:
                    self.pair_selection_update()

        if self.animation:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower') # , vmin=-1, vmax=1)
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

    def run_tests(self, n_tests=20, steps_per_measurement=10, measurments_per_test=50):  # 20 10 100
        T_arr = np.linspace(0.1, 3, n_tests)
        self.c = 0.8  # Fixed value given in the question.
        # Arrays to store values for plots
        heat_capacities = np.zeros(n_tests)
        errors = np.zeros(n_tests)
        suceptabilities = np.zeros(n_tests)
        magnetisations = np.zeros(n_tests)
        energies = np.zeros(n_tests)
        for i in range(n_tests):
            self.T = T_arr[i]
            self.grid = np.random.choice([1, -1, 0], size=(self.N,self.N), p=[self.c*0.5,self.c*0.5,1-self.c])
            sys.stdout.write("Simulation progress: %.1f%%\r"
                            % ((100 * i / n_tests)))
            sys.stdout.flush()  # Prints progress of simulation.

            # Stabalising sweeps
            self.sweeps_per_update = 100  # 100
            self.update(1)
            # Reset the sweeps per update.
            self.sweeps_per_update = steps_per_measurement

            E_arr = np.zeros(measurments_per_test)
            M_arr = np.zeros(measurments_per_test)
            for j in range(measurments_per_test):
                E_arr[j] = self.sys_energy()
                M_arr[j] = self.magnetisation()
                self.update(1)
            heat_capacities[i] = self.heat_cap(E_arr)
            errors[i] = self.bootstrap_errors(E_arr)
            suceptabilities[i] = self.suceptability(M_arr)
            magnetisations[i] = np.average(M_arr)
            energies[i] = np.average(E_arr)
        heat_capacities = np.array(list(zip(T_arr, heat_capacities, errors)))
        suceptabilities = np.array(list(zip(T_arr, suceptabilities)))
        magnetisations = np.array(list(zip(T_arr, magnetisations)))
        energies = np.array(list(zip(T_arr, energies)))

        np.savetxt('heat_capacities', heat_capacities)
        np.savetxt('suceptabilities', suceptabilities)
        np.savetxt('magnetisations', magnetisations)
        np.savetxt('energies', energies)

    def magnetisation(self):
        return np.abs(np.sum(self.grid))

    def sys_energy(self):
        E = 0
        for i in range(self.N):
            for j in range(self.N):
                # Only consider N and E components to stop double counting.
                E += self.J*(self.grid[(i-1+self.N)%self.N][j] + self.grid[i][(j-1+self.N)%self.N]) * self.grid[i][j]  # Removed a -ve
        return E

    def heat_cap(self, E_arr):
        norm_factor = 1/(self.kb*self.T**2)
        Cv = np.average(np.square(E_arr)) - np.square(np.average(E_arr))
        Cv *= norm_factor
        return Cv

    def suceptability(self, M_arr):
        print(M_arr)
        norm_factor = 1/(self.kb*self.T**2)
        chi = np.average(np.square(M_arr)) - np.square(np.average(M_arr))
        chi *= norm_factor
        return chi

    def bootstrap_errors(self, E_arr, k=100):
        heat_cap = np.zeros(k)

        for i in range(k):
            rand_selection = np.random.choice(E_arr, len(E_arr))
            heat_cap[i] = self.heat_cap(rand_selection)
        sigma = np.sqrt(np.average(np.square(heat_cap))
                 - np.square(np.average(heat_cap)))
        return sigma
