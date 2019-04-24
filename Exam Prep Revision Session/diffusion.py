import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class Diffusion:
    def __init__(self, N, mean_val=0.5, sigma=10, kapa=0.01, D=1, dx=1, dt=0.01):
        self.N = N
        self.sigma = sigma
        self.kapa = kapa
        self.D = D
        self.dx = dx
        self.dt = dt
        self.sweeps_per_update = 1

        noise = 0.1
        self.psi_grid = np.random.uniform(low=mean_val-noise, high=mean_val+noise, size=(N, N))
        self.next_psi_grid = np.zeros((N, N))
        self.rho_grid = np.zeros((N, N))
        self.init_pho_grid()
        self.animation = True

    def init_pho_grid(self):
        """
        Sets up the rho grid such that it follows the distribution,
            rho = exp(-r^2/sigma^2)
        about the centre of the array.
        """
        centre_x = int(self.N/2)
        centre_y = int(self.N/2)
        for i in range(self.N):
            for j in range(self.N):
                r = np.sqrt((i - centre_y)**2 + (j - centre_x)**2)
                self.rho_grid[i][j] = np.exp(-r**2/ self.sigma**2)
        # Quick check to make sure that this works.
        # plt.imshow(self.rho_grid)
        # plt.show()

    def laplacian(self, i, j):
        l1 = self.psi_grid[(i+1+self.N) % self.N][j] + self.psi_grid[(i-1+self.N) % self.N][j]
        l2 = self.psi_grid[i][(j+1+self.N) % self.N] + self.psi_grid[i][(j-1+self.N) % self.N]
        l3 = -4*self.psi_grid[i][j]
        return l1 + l2 + l3

    def rolling_nn_sum(self):
        self.rolling_nn = (np.roll(self.psi_grid, 1, 0) + np.roll(self.psi_grid, -1, 0)
            + np.roll(self.psi_grid, 0, 1) + np.roll(self.psi_grid, 0, -1))

    def rolling_laplacian(self, i, j):
        return self.rolling_nn[i][j] - 4*self.psi_grid[i][j]

    def update(self, k):
        for z in range(self.sweeps_per_update):
            # rand_choices = np.random.randint(low=0, high=self.N, size=(self.N, self.N, 2))
            # update_probabilities = np.random.uniform(size=(self.N, self.N))
            self.rolling_nn_sum()
            for i in range(self.N):
                for j in range(self.N):
                    # self.next_psi_grid[i][j] = (self.psi_grid[i][j] + self.D * self.dt / self.dx**2
                    #     * self.laplacian(i, j) + self.rho_grid[i][j] - self.kapa * self.psi_grid[i][j])
                    self.next_psi_grid[i][j] = (self.psi_grid[i][j] + self.D * self.dt / self.dx**2
                        * self.rolling_laplacian(i, j) + self.rho_grid[i][j] - self.kapa * self.psi_grid[i][j])
        self.psi_grid = self.next_psi_grid.copy()
        if self.animation:
            self.fig.clear()
            plt.imshow(self.psi_grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower')# , vmin=0, vmax=2)
            plt.colorbar()

    def measure_mean_psi(self, steps=5000):
        mean_psi = np.zeros(steps)
        self.sweeps_per_update = 1
        for i in range(steps):
            mean_psi[i] = np.mean(self.psi_grid)
            self.update(1)
        lin = [x for x in range(steps)]
        mean_psi = np.array(list(zip(mean_psi, lin)))
        np.savetxt("Mean_psi_val_with_time.txt", mean_psi, header='The mean psi value changing with time over %d steps.' % steps)

    def psi_with_r(self):
        # 5000 form the leveling off of the graph from part 3.
        self.sweeps_per_update = 5000
        self.update(1)
        centre_x = int(self.N/2)
        centre_y = int(self.N/2)
        # results = np.zeros(self.N**2, 2)
        r_vals = np.zeros(self.N**2)
        psi_vals = np.zeros(self.N**2)
        z = 0
        for i in range(self.N):
            for j in range(self.N):
                r_vals[z] = np.sqrt((i - centre_y)**2 + (j - centre_x)**2)
                psi_vals[z] = self.psi_grid[i][j]
                z += 1
        results = np.array(list(zip(psi_vals, r_vals)))
        np.savetxt("psi_with_r.txt", results, header="psi varying with r after %d equilibrium sweeps." % self.sweeps_per_update)

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
        anim = FuncAnimation(self.fig, self.update, interval=5)
        plt.show()
