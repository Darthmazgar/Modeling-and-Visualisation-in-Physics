import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import exp
import sys

class Other: ################
    def __init__(self, N, J=-1, kb=1, T=1):
        self.N = N
        self.grid = np.random.choice([-1, 1], size=(N,N))
        self.animation = False
        self.sweeps_per_update = 1

    def update(self, k):
        for z in range(self.sweeps_per_update):
            pass
        if self.animation:
            self.fig.clear()
            plt.imshow(self.grid, interpolation='nearest',
                           cmap='coolwarm', origin='lower' , vmin=-1, vmax=1)
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
