import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class Sirs:
    def __init__(self, N, M, p1, p2, p3, anim=True):
        self.N = N
        self.M = M
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        # 0: suceptable, 1: infected, 2: recovered.
        self.grid = np.random.choice([0, 1], size=(N, M), p=[1/2, 1/2])
        self.sweeps = 1
        if anim:
            self.fig = plt.figure()

    def check_nn(self, i, j):
        if self.grid[(i + 1 +self.N) % self.N][j] == 1:
            return 1
        elif self.grid[(i - 1 +self.N) % self.N][j] == 1:
            return 1
        elif self.grid[i][(j + 1 + self.M) % self.M] == 1:
            return 1
        elif self.grid[i][(j - 1 + self.M) % self.M] == 1:
            return 1
        else:
            return 0

    def update(self, k, anim=True):
        for k in range(self.sweeps):
            rand_ar = np.random.random((self.N,self.M))
            for i in range(self.N):
                for j in range(self.M):
                    state = self.grid[i][j]
                    if state == 0:
                        p = self.p1
                    elif state == 1:
                        p = self.p2
                    elif state == 2:
                        p = self.p3
                    if state == 0:
                        if self.check_nn(i, j) and rand_ar[i][j] <= p:
                            self.grid[i][j] = (self.grid[i][j] + 1) % 3
                    else:
                        if rand_ar[i][j] <= p:
                            self.grid[i][j] = (self.grid[i][j] + 1) % 3
        if anim:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                       cmap='Blues', vmin=0, vmax=2)
        return self.grid

    def run_animation(self):
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


def main(argv):
    if len(argv) != 5:
        print("gol.py M N P1 P2 P3")
        sys.exit()
    M = int(argv[0])
    N = int(argv[1])
    P1 = float(argv[2])
    P2 = float(argv[3])
    P3 = float(argv[4])

    s = Sirs(M, N, P1, P2, P3)
    s.run_animation()

main(sys.argv[1:])
