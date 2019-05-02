import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import exp
import sys

class Ising:
    def __init__(self, N, J=-1, kb=1, T=1):
        self.N = N
        self.J = J
        self.kb = kb
        self.T = T
        self.h = 100
        self.grid = np.random.choice([-1, 1], size=(N,N))
        self.sum_nearest_neighbours = None
        self.animation = False
        self.sweeps_per_update = 10

    def rolling_nn_sum(self):
        self.sum_nearest_neighbours = np.roll(self.grid, 1,0) + np.roll(self.grid, -1, 0) + np.roll(self.grid, 1, 1) + np.roll(self.grid, -1, 1)
        return self.sum_nearest_neighbours

    def dE(self, i, j):
        Ei = -self.J*self.sum_nearest_neighbours[i][j]*self.grid[i][j] - self.h* self.grid[i][j]
        Ef = -self.J*self.sum_nearest_neighbours[i][j]*(-1*self.grid[i][j]) - self.h*(-1*self.grid[i][j])
        dE = Ef - Ei
        return dE

    def P(self, dE):
        p = np.exp(-dE / (self.kb * self.T))
        # print(p)
        return p

    def magnetisation(self):
        m = np.sum(self.grid)
        return m

    def staggered_mag(self):
        sum = 0
        for i in range(self.N):
            for j in range(self.N):
                sgn = -1**(i+j)
                sum += sgn * self.grid[i][j]
        return sum

    def sys_energy(self):
        eng = 0
        for i in range(self.N):
            for j in range(self.N):
                eng += self.E(i, j) / 2  # /2 to not have double counting of nearest neighbours.
        return eng

    def E(self, i, j):
        E = -self.J*self.sum_nearest_neighbours[i][j]*self.grid[i][j] - self.h* self.grid[i][j]
        return E

    def set_h(self, x, y, n):
        P = 25
        tau = 10000
        h0 = 10
        self.h = h0 * np.cos(2*pi*x/P) * np.cos(2*np.pi*y/P)*np.sin(2*np.pi*n/tau)

    def update(self, k):
        for z in range(self.sweeps_per_update):
            self.rolling_nn_sum()
            rand_no = np.random.uniform(size=(self.N, self.N))
            update_positions = np.random.randint(self.N, size=(self.N, self.N, 2))
            for i in range(self.N):
                for j in range(self.N):
                    # E = self.E(i, j)
                    dE = self.dE(update_positions[i][j][0], update_positions[i][j][1])
                    if dE <= 0 or self.P(dE) < rand_no[i][j]:
                        # self.grid[i][j] *= -1
                        self.grid[update_positions[i][j][0]][update_positions[i][j][1]] *= -1
                    else:
                        continue
        if self.animation:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower' , vmin=-1, vmax=1)
            plt.colorbar()

    def run_tests(self, n_steps=20, steps_per_test=10, no_tests=500):  # 20, 10, 100
        hs = np.linspace(0, 10, n_steps)
        mag_mean = np.zeros(n_steps)
        mag_var = np.zeros(n_steps)
        stag_mag_mean = np.zeros(n_steps)
        stag_mag_var = np.zeros(n_steps)
        mean_eng = np.zeros(n_steps)

        for i in range(n_steps):
            sys.stdout.write("Simulation progress: %.1f%%\r"
                            % ((100 * i / n_steps)))
            sys.stdout.flush()  # Prints progress of simulation.

            # Reset the grid between h changes.
            self.grid = np.random.choice([-1, 1], size=(self.N,self.N))

            mag = np.zeros(no_tests)
            stag_mag = np.zeros(no_tests)
            eng = np.zeros(no_tests)

            # Move on to the next h value
            self.h = hs[i]

            self.sweeps_per_update = 500  # Equilibriation sweeps  # 500
            self.update(1)  # Run equilibriation

            self.sweeps_per_update = steps_per_test

            for j in range(no_tests):
                # Make measurements
                mag[j] = self.magnetisation()
                stag_mag[j] = self.staggered_mag()
                eng[j] = self.sys_energy()

                # Update steps_per_test times
                self.update(1)

            # Calculate averages and variences
            mag_mean[i] = np.mean(mag)
            mag_var[i] = np.var(mag)
            stag_mag_mean[i] = np.mean(stag_mag)
            stag_mag_var[i] = np.var(stag_mag)
            mean_eng[i] = np.mean(eng)

        # Output files
        mag_mean_zipped = np.array(list(zip(hs, mag_mean)))
        mag_var_zipped = np.array(list(zip(hs, mag_var)))
        stag_mag_mean_zipped = np.array(list(zip(hs, stag_mag_mean)))
        stag_mag_var_zipped = np.array(list(zip(hs, stag_mag_var)))
        mean_eng_zipped = np.array(list(zip(hs, mean_eng)))

        np.savetxt('mag_mean_zipped', mag_mean_zipped)
        np.savetxt('mag_var_zipped', mag_var_zipped)
        np.savetxt('stag_mag_mean_zipped', stag_mag_mean_zipped)
        np.savetxt('stag_mag_var_zipped', stag_mag_var_zipped)
        np.savetxt('mean_eng_zipped', mean_eng_zipped)


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
