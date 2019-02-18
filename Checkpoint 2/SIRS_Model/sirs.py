import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

class Sirs:
    def __init__(self, N, M, p1, p2, p3, immune=0, test=False, anim=True):
        self.N = N
        self.M = M
        self.p1 = p1
        self.p2 = p2
        if test:
            self.p2 = 0.5
        self.p3 = p3
        self.test = test
        # 0: suceptable, 1: infected, 2: recovered.
        frac = (1 - immune) / 2.
        self.grid = np.random.choice([0, 1, 3], size=(N, M),
                    p=[frac, frac, immune])
        if anim:
            self.fig = plt.figure()

    def check_nn(self, i, j):
        if self.grid[(i + 1 +self.N) % self.N][j] == 1:
            return True
        elif self.grid[(i - 1 +self.N) % self.N][j] == 1:
            return True
        elif self.grid[i][(j + 1 + self.M) % self.M] == 1:
            return True
        elif self.grid[i][(j - 1 + self.M) % self.M] == 1:
            return True
        else:
            return False

    def update(self, k, sweeps=1, anim=True):
        for k in range(sweeps):
            rand_ar = np.random.random((self.N,self.M))  # Gen probabilities.
            for i in range(self.N):
                for j in range(self.M):
                    state = self.grid[i][j]
                    if state == 3:  # If immune then leave
                        continue
                    elif state == 0:  # If suceptable
                        p = self.p1
                    elif state == 1:  # If Infected
                        p = self.p2
                    elif state == 2:  # If recovering
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
                       cmap='brg', vmin=0, vmax=2)
        return self.grid

    def measure_infected(self):
        infected = 0
        for i in range(self.N):
            for j in range(self.M):
                if self.grid[i][j] == 1:
                    infected += 1
        return infected / (self.N * self.M)

    def run_test(self, resolution=10, sweeps_per_test=100,
                measurements_per_test=10, show=True, save=True):
        heatmap = np.zeros((resolution, resolution))  # , measurements_per_test))
        p1_ar = np.linspace(0, 1, resolution)
        p3_ar = np.linspace(0, 1, resolution)
        self.sweeps_per_test = sweeps_per_test
        for i in range(resolution):
            self.p1 = p1_ar[i]
            for j in range(resolution):
                sys.stdout.write("Simulation progress: %.1f%%: Part progress: %.1f%%\r"
                                % ((100 * i / resolution), (100 * j / resolution)))
                sys.stdout.flush()  # Prints progress of simulation.
                self.p3 = p3_ar[j]
                test_results = np.zeros(measurements_per_test)
                for k in range(measurements_per_test):
                    self.grid = np.random.choice([0, 1], size=(self.N, self.M),
                                                p=[1/2, 1/2])
                    self.update(1, sweeps=sweeps_per_test, anim=False)
                    test_results[k] = self.measure_infected()
                    # heatmap[i][j][k] = test_results[k]
                heatmap[i][j] = np.average(test_results)
        if show:
            plt.imshow(heatmap, interpolation='nearest',
                       cmap='brg', vmin=0, vmax=2)
            plt.xlabel("P3")
            plt.ylabel("P1")
            plt.title("Sirs test")
            plt.show()
        if save:
            np.savetxt('sirs_heatmap.txt', heatmap)

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

        anim = FuncAnimation(self.fig, self.update, interval=5, frames=50)
        plt.show()


def main(argv):
    if len(argv) != 7:
        print("gol.py M N P1 P2 P3 Pimmune anim:0 / test:1")
        sys.exit()
    M = int(argv[0])
    N = int(argv[1])
    P1 = float(argv[2])
    P2 = float(argv[3])
    P3 = float(argv[4])
    P4 = float(argv[5])

    if argv[6] == '0' or argv[6] == 'anim':
        s = Sirs(M, N, P1, P2, P3, immune=P4)
        s.run_animation()
    elif argv[6] == '1' or argv[6] == 'test':
        s = Sirs(M, N, P1, P2, P3, immune=P4, test=True)
        s.run_test()

main(sys.argv[1:])
