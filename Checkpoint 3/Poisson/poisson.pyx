import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys


class Poisson:
    def __init__(self):
        pass

    def update(self, int k):
        pass

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




if __name__ == '__main__':
    main(sys.argv[1:])
