import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys

class CahnHil:
    def __init__(self, N, dx=.7, dt=.7, m=.1, k=.1, a=.1, b=.1, phi_0=0,
                anim=True, f_test=False):
        print("Running pyx")
        self.N = N
        self.dt = dt
        self.dx = dx
        self.M = m
        self.k = k
        self.a = a
        self.b = b
        self.steps_per_update = 50
        flucs = 0.01
        self.grid = np.random.uniform(low=-flucs+phi_0, high=flucs+phi_0,
                                    size=(self.N, self.N))
        # self.grid = np.ones((self.N, self.N)) * -1
        self.next_grid = np.zeros((self.N, self.N))
        # self.add_drop(5, int(N/2), int(N/2))
        if anim:
            self.fig = plt.figure()

    def update(self, int k, anim=True):
        cdef int z, i, j
        if k != 0 and anim != 'run':
            print("\nUpdate: %d" % k)
            print("Frames: %d" % (k * self.steps_per_update))

        for z in range(self.steps_per_update):
            if anim == 'run':
                sys.stdout.write("Simulation progress: %.1f%%\r"
                                % ((100 * z / self.steps_per_update)))
                sys.stdout.flush()  # Prints progress of simulation.
            for i in range(self.N):
                for j in range(self.N):
                    self.next_grid[i][j] = self.next_phi(i, j)
            self.grid = self.next_grid.copy()
        if anim:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                       cmap='coolwarm', vmin=-1, vmax=1, origin='lower')
        return self.grid

    def display(self, sweeps):
        plt.imshow(self.grid, interpolation='nearest',
                   cmap='coolwarm', vmin=-1, vmax=1, origin='lower')
        plt.show()

    def add_drop(self, int size, int x, int y):
        cdef int i, j, dens, r
        r = int(size/2)
        for i in range(x-r, x+r):
            for j in range(y-r, y+r):
                # dens = (np.abs(r-i) + np.abs(r-j)) / 2*r
                dens=1
                self.grid[(i+ self.N) % self.N][(j+self.N) % self.N] = dens

    def mu(self, int i, int j):
        cdef double t1, norm_fact, t2
        t1 = (-self.a * self.grid[(i + self.N) % self.N][(j + self.N) % self.N]
            + self.b * self.grid[(i + self.N) % self.N][(j + self.N) % self.N]**3)
        norm_factor = - (self.k / self.dx**2)
        t2 = (self.grid[(i+1 + self.N) % self.N][(j + self.N) % self.N]
            + self.grid[(i-1 + self.N) % self.N][(j + self.N) % self.N]
            + self.grid[(i + self.N) % self.N][(j+1 + self.N) % self.N]
            + self.grid[(i + self.N) % self.N][(j-1) + self.N % self.N]
            - 4 * self.grid[(i + self.N) % self.N][(j + self.N) % self.N])
        return t1 + norm_factor * t2

    def next_phi(self, int i, int j):
        cdef double norm_fact, expression, new_phi, next_val
        norm_factor = (self.M * self.dt) / self.dx**2
        expression = (self.mu(i+1, j) + self.mu(i-1, j) + self.mu(i, j+1)
                    + self.mu(i, j-1) - 4 * self.mu(i, j))
        new_phi = norm_factor * expression
        next_val = self.grid[i][j] + new_phi
        # print(self.grid[i][j])
        return next_val

    def free_eng(self, int i, int j):
        cdef double t1, t2, t3, f
        t1 = - (self.a/2) * self.grid[i][j]**2
        t2 = (self.a/4) * self.grid[i][j]**4
        t3 = (self.k/2) * (self.grid[(i+1 + self.N) % self.N][j]
            + self.grid[(i-1 + self.N) % self.N][j]
            + self.grid[i][(j+1 + self.N) % self.N]
            + self.grid[i][(j-1) + self.N % self.N]
            - 4 * self.grid[i][j])**2
        f = t1 + t2 + t3
        return f

    def run_f_test(self, int test_length, save=False):
        cdef int z, i, j
        f_engs = np.zeros(test_length)
        for z in range(test_length):
            sys.stdout.write("Simulation progress: %.1f%%\r"
                            % ((100 * z / test_length)))
            sys.stdout.flush()  # Prints progress of simulation.
            for i in range(self.N):
                for j in range(self.N):
                    f_engs[z] += self.free_eng(i, j)
            f_engs[z] /= self.N**2
            self.update(1, anim=False)
            print(f_engs[z])
        np.savetxt('free_energy.txt', f_engs, header='Free energy with time'
            +' for a time step of %.2f, space step %.2f' % (self.dt, self.dx))
        xs = np.linspace(0, test_length*self.dt*self.steps_per_update,
                        test_length)
        plt.plot(xs, f_engs)
        plt.xlabel('Time')
        plt.ylabel("Free energy")
        plt.title("Free energy with time using the Cahn Hilliard equation")
        if save:
            plt.savefig('Free-energy-with-time.png')
        plt.show()

    def animate(self):
        anim = FuncAnimation(self.fig, self.update)
        plt.show()

    def run_animation(self):
        """
        Gives the ability to click on the animation canvas to play and pause.
        """
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
        anim = FuncAnimation(self.fig, self.update, interval=5)
        plt.show()


def main(args):
    if len(args) != 5:
        print("\nTo few arguments.\n")
        print("python cahn_hilliard.py N dt dx phi_0 0:anim or 1:test")
        sys.exit()
    N = int(args[0])
    dt = float(args[1])
    dx = float(args[2])
    phi_0 = float(args[3])
    test = args[4]

    if args[4] == '0' or args[4] == 'anim':
        ch = CahnHil(N, dt=dt, dx=dx, phi_0=phi_0)
        if len(args) == 6:
            for i in range(int(args[5])):
                ch.add_drop(np.random.uniform()*15, np.random.randint(N), np.random.randint(N))
        ch.run_animation()
    elif args[4] == '1' or args[4] == 'test':
        cd = CahnHil(N, dt=dt, dx=dx, f_test=True,phi_0=phi_0)
        cd.run_f_test(1000)


if __name__ == '__main__':
    main(sys.argv[1:])
