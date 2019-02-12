import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class GameOfLife:
    def __init__(self, N, M, dens, set_state='rand', anim=True):
        self.N = N
        self.M = M
        self.sweeps = 1
        self.set_state = set_state
        if set_state == 'rand':
            self.grid = np.random.choice([0, 1], size=(N, M), p=[1 - dens, dens])
        elif set_state == 'oscilator':
            self.grid = np.zeros((N, M))
            self.make_oscilator(int(N/2), int(M/2))
        elif set_state == 'glider':
            self.grid = np.zeros((N, M))
            self.make_glider(int(N/2), int(M/2))
        elif set_state == 'glider gun':
            self.grid = np.zeros((N, M))
            self.make_glider_gun(50, 50)
        self.new_grid = self.grid.copy()  # Have to use .copy!
        if anim:
            self.fig = plt.figure()
        # self.init_kaw_grid()

    def init_half_grid(self):
        """
        Initialise the grid with one half all 1 and the other all 0.
        """
        ones = np.ones((int(self.N / 2), self.M))
        neg_ones = ones * -1
        self.grid = np.concatenate((ones, neg_ones))
        return self.grid

    def make_oscilator(self, x, y):
        self.grid[(x + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + 1 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x - 1 + self.M) % self.M][(y + self.M) % self.M] = 1

    def make_glider(self, x, y):
        self.grid[(x + 1 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        self.grid[(x + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 1 + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        self.grid[(x - 1 + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        return 1

    def make_glider_gun(self, x, y):
        # Left square
        self.grid[(x + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + 1 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 1 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1

        # Right square
        self.grid[(x + 34 + self.N) % self.N][(y + 2 + self.M) % self.M] = 1
        self.grid[(x + 34 + 1 + self.N) % self.N][(y + 2 + self.M) % self.M] = 1
        self.grid[(x + 34 + self.N) % self.N][(y + 2 - 1 + self.M) % self.M] = 1
        self.grid[(x + 34 + 1 + self.N) % self.N][(y + 2 - 1 + self.M) % self.M] = 1

        # Left bit
        self.grid[(x + 10 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + 10 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 10 + self.N) % self.N][(y - 2 + self.M) % self.M] = 1
        self.grid[(x + 11 + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        self.grid[(x + 11 + self.N) % self.N][(y - 3 + self.M) % self.M] = 1
        self.grid[(x + 12 + self.N) % self.N][(y + 2 + self.M) % self.M] = 1
        self.grid[(x + 12 + self.N) % self.N][(y - 4 + self.M) % self.M] = 1
        self.grid[(x + 13 + self.N) % self.N][(y + 2 + self.M) % self.M] = 1
        self.grid[(x + 13 + self.N) % self.N][(y - 4 + self.M) % self.M] = 1
        self.grid[(x + 14 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 15 + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        self.grid[(x + 15 + self.N) % self.N][(y - 3 + self.M) % self.M] = 1
        self.grid[(x + 16 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + 16 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 16 + self.N) % self.N][(y - 2 + self.M) % self.M] = 1
        self.grid[(x + 17 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1

        # Right bit
        self.grid[(x + 20 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + 20 + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        self.grid[(x + 20 + self.N) % self.N][(y + 2 + self.M) % self.M] = 1
        self.grid[(x + 21 + self.N) % self.N][(y + self.M) % self.M] = 1
        self.grid[(x + 21 + self.N) % self.N][(y + 1 + self.M) % self.M] = 1
        self.grid[(x + 21 + self.N) % self.N][(y + 2 + self.M) % self.M] = 1
        self.grid[(x + 22 + self.N) % self.N][(y + 3 + self.M) % self.M] = 1
        self.grid[(x + 22 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 24 + self.N) % self.N][(y + 3 + self.M) % self.M] = 1
        self.grid[(x + 24 + self.N) % self.N][(y + 4 + self.M) % self.M] = 1
        self.grid[(x + 24 + self.N) % self.N][(y - 1 + self.M) % self.M] = 1
        self.grid[(x + 24 + self.N) % self.N][(y - 2 + self.M) % self.M] = 1
        return 1

    def com_tracking(self, s_x=0, s_y=0):
        s_x=int(self.N/2)
        s_y=int(self.M/2)
        n = np.sum(self.grid)  # count alive cells so it dosent just have to be a glider.
        sum = 0
        for i in range(self.N):
            for j in range(self.M):
                state = self.grid[i][j]
                if state:
                    sum += np.sqrt((i - s_x)**2 + (j - s_y)**2)
        r = sum / n
        return r

    def update(self, k, anim=True):
        for z in range(self.sweeps):
            for i in range(self.N):
                for j in range(self.M):
                    count = 0
                    state = self.grid[i][j]
                    for x in range(-1 + i, 2 + i):
                        for y in range(-1 + j, 2 + j):
                            if x == i and y == j:
                                continue
                            count += self.grid[(x + self.N) % self.N]\
                                              [(y + self.M) % self.M]
                    if state == 0 and count == 3:
                        # print("Born: %d" % count)
                        self.new_grid[i][j] = 1
                    elif state == 1 and (count < 2 or count > 3):
                        # print("Died: %d" % count)
                        self.new_grid[i][j] = 0
            self.grid = self.new_grid.copy()  # Have to use .copy!
        if self.set_state == 'glider':
            r = self.com_tracking()
            print(self.com_tracking())
        if anim:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='None',
                       cmap='Blues', vmin=0, vmax=1)


        return self.grid

    def population(self, pr=True):
        size = self.N * self.M
        pop = np.sum(self.grid)
        frac = pop / size
        if pr:
            print("Alive population: %.1f%%" % frac*100)
        return frac

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
                anim.event_source.stop()
                anim_running = False
            else:
                anim.event_source.start()
                anim_running = True

        self.fig.canvas.mpl_connect('button_press_event', onClick)

        anim = FuncAnimation(self.fig, self.update)
        plt.show()


def main(argv):
    if len(argv) != 5:
        print("gol.py M N Density(0->1) Setup(rand:0, glider:1, oscilator:2)\
              Test(anim:0, test:1)")
        sys.exit()
    M = int(argv[0])
    N = int(argv[1])
    dens = float(argv[2])

    if argv[3] == '0' or argv[3] == 'rand':
        gol = GameOfLife(N, M, dens)
    elif argv[3] == '1' or argv[3] == 'glider':
        gol = GameOfLife(N, M, dens, set_state='glider')
    elif argv[3] == '2' or argv[3] == 'oscilator':
        gol = GameOfLife(N, M, dens, set_state='oscilator')
    elif argv[3] == '3' or argv[3] == 'glider_gun':
        gol = GameOfLife(N, M, dens, set_state='glider gun')
    if argv[4] == '0' or argv[4] == 'anim':
        # gol.animate()
        gol.run_animation()
    elif argv[4] == '1' or argv[4] == 'test':
        gol = GameOfLife(N, M, dens, anim=False)
    else:
        print("Not valid input for anim or test.")



main(sys.argv[1:])
