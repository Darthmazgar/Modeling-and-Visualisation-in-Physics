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
    s = Sirs(5, 5, 0.1, 0.1, 0.1)

main(sys.argv[1:])
