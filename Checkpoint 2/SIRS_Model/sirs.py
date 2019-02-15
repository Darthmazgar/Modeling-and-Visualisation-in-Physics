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
        self.grid = np.random.choice([0, 1, 2], size=(N, M), p=[1/3, 1/3, 1/3])
        if anim:
            self.fig = plt.figure()

    def update(self, k, anim=True):

        if anim:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='None',
                       cmap='Blues', vmin=0, vmax=1)
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

    s = Sirs(5, 5, 0.1, 0.1, 0.1)

main(sys.argv[1:])
