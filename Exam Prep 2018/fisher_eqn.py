import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class FisherEqn:
    def __init__(self, N, R, dt, dx=1, D=1, alpha=1):
        self.N = N
        self.R = R
        self.D = D
        self.dx = dx
        self.dt = dt
        self.alpha = alpha
        self.phi_grid = np.zeros((N, N))
        self.init_grid()
        self.sweeps_per_update = 100

        self.animation = True

    def init_grid(self):
        """
        Adds values of 1 to the phi grid for all points in a radius |r| < R.
        Due to boundary conds this is added at the middle of the grid.
        """
        ox = int(self.N/2)
        oy = int(self.N/2)
        # Only working over a small square which encloses all points where |r| < R.
        for i in range(-int(self.R) + oy + 1, int(self.R) + oy + 1):
            for j in range(-int(self.R) + ox + 1, int(self.R) + ox + 1):
                # Check if |r| < R.
                mag_r = np.sqrt((j - ox)**2 + (i - oy)**2)
                if(mag_r < self.R):
                    self.phi_grid[i][j] = 1

    def del_sq(self, i, j):
        """ Applies the laplacian operator using periodic boundary conditions
        on the point (i, j) in the phi_grid.
        """
        l  = self.dt / self.dx**2
        l1 = self.phi_grid[(i + 1 + self.N) % self.N][j] + self.phi_grid[(i - 1 + self.N) % self.N][j]
        l2 = self.phi_grid[i][(j + 1 + self.N) % self.N] + self.phi_grid[i][(j - 1 + self.N) % self.N]
        l3 = - 4 * self.phi_grid[i][j]
        return self.phi_grid[i][j] + l * (l1 + l2 +l3)


    def update(self, k):
        """
        Updates the phi_grid according to the Fisher Equation.
        d psi / dt = D * Del^2 psi + alpha * psi(1 - psi)
        """
        for z in range(self.sweeps_per_update):
            next_grid = np.zeros((self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    p1 = self.D * self.del_sq(i, j) * self.phi_grid[i][j]
                    p2 = self.alpha * self.phi_grid[i][j] * (1 - self.phi_grid[i][j])
                    next_grid[i][j] = p1 + p2

            self.phi_grid = next_grid.copy()

        if self.animation:
            self.fig.clear()
            # plt.title("Potential Strength")
            plt.imshow(self.phi_grid, interpolation='nearest',
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
