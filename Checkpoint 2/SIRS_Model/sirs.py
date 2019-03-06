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
            rand_xs = np.random.randint(self.N, size=self.N)
            rand_ys = np.random.randint(self.M, size=self.M)
            for i in range(self.N):
                for j in range(self.M):
                    x = rand_xs[i]
                    y = rand_ys[j]
                    # x=i
                    # y=j
                    state = self.grid[x][y]
                    if state == 3:  # If immune then leave
                        continue
                    elif state == 0:  # If suceptable
                        p = self.p1
                    elif state == 1:  # If Infected
                        p = self.p2
                    elif state == 2:  # If recovering
                        p = self.p3
                    if state == 0:
                        if self.check_nn(x, y) and rand_ar[i][j] <= p:
                            self.grid[x][y] = (self.grid[x][y] + 1) % 3
                    else:
                        if rand_ar[x][y] <= p:
                            self.grid[x][y] = (self.grid[x][y] + 1) % 3
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

    def phase_test(self, resolution=20, sweeps_per_test=25,
                measurements_per_test=2, show=True, save=True):
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
            plt.imshow(heatmap, interpolation='nearest', cmap='brg', extent=[0,1,0,1], origin='lower')
            plt.xlabel("P3")
            plt.ylabel("P1")
            plt.title("Sirs test")
            plt.show()
        if save:
            np.savetxt('sirs_heatmap1.txt', heatmap)

    def contour_test(self, resolution=10, sweeps_per_test=50,
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
                self.grid = np.random.choice([0, 1], size=(self.N, self.M),
                                            p=[1/2, 1/2])
                self.update(1, sweeps=100, anim=False)
                for k in range(measurements_per_test):
                    self.update(1, sweeps=sweeps_per_test, anim=False)
                    test_results[k] = self.measure_infected()
                heatmap[i][j] = (np.average(np.square(test_results)) - np.square(np.average(test_results))) / (self.N * self.M)
        if show:
            plt.imshow(heatmap, interpolation='nearest', cmap='brg', extent=[0,1,0,1], origin='lower') # , vmin=0,
            plt.xlabel("P1")
            plt.ylabel("P3")
            plt.title("Contour test")
            plt.show()
        if save:
            np.savetxt('contour_test.txt', heatmap)

    def slice_test(self, resolution=50, sweeps_per_test=10,
                measurements_per_test=200, show=True, save=True):
        p1_ar = np.linspace(0.2, 0.5, resolution)
        var = np.zeros(resolution)
        self.sweeps_per_test = sweeps_per_test
        for i in range(resolution):
            self.p1 = p1_ar[i]  # Set p1 value to next.

            # sys.stdout.write("Simulation progress: %.1f%%\r"
            #                 % ((100 * i / resolution)))
            # sys.stdout.flush()  # Prints progress of simulation.

            test_results = np.zeros(measurements_per_test)
            self.grid = np.random.choice([0, 1], size=(self.N, self.M),  # Reset the grid
                                        p=[1/2, 1/2])
            self.update(1, sweeps=100, anim=False)  # Let sys settle
            for k in range(measurements_per_test):
                sys.stdout.write("Simulation progress: %.1f%%: Part progress: %.1f%%\r"
                                % ((100 * i / resolution), (100 * k / measurements_per_test)))
                sys.stdout.flush()  # Prints progress of simulation.
                self.update(1, sweeps=sweeps_per_test, anim=False)  # Run for so many sweeps
                test_results[k] = self.measure_infected()  # measure infected
            var[i] = (np.average(np.square(test_results)) - np.square(np.average(test_results))) / (self.N * self.M)
        if show:
            plt.plot(p1_ar, var)
            plt.xlabel("P1")
            plt.ylabel("Varience")
            plt.title("Slice test (P2=0.5, P3=0.5)")
            plt.show()
        if save:
            np.savetxt('contour_test.txt', var)

    def waves_test(self, resolution=250, sweeps_per_test=10,
                    show=True, save=True):
        infected_ar = np.zeros(resolution)
        sweeps_ar = np.linspace(0, resolution*sweeps_per_test, resolution)
        self.update(1, sweeps=10, anim=False)
        for i in range(resolution):
            sys.stdout.write("Simulation progress: %.1f%%\r"
                            % ((100 * i / resolution)))
            sys.stdout.flush()  # Prints progress of simulation.
            self.update(1, sweeps=sweeps_per_test, anim=False)
            infected_ar[i] = self.measure_infected()
        print(infected_ar)
        plt.plot(sweeps_ar, infected_ar)
        plt.xlabel("Sweeps")
        plt.ylabel("Infected fraction")
        plt.title("Infected fraction with time")
        if save:
            np.savetxt("infected_frac_time.txt", infected_ar)
        if show:
            plt.show()

    def immune_test(self, resolution=100, sweeps_per_test=1000, measurements_per_test=2,  # res: 50, spt: 100, mpt: 10
                    show=True, save=True):
        immune_ar = np.linspace(0, 1, resolution)
        ys = np.zeros(resolution)
        yerr = np.zeros(resolution)
        for i in range(resolution):
            test_results = np.zeros(measurements_per_test)
            sys.stdout.write("Simulation progress: %.1f%%\r" % ((100 * i / resolution)))
            sys.stdout.flush()  # Prints progress of simulation.

            for j in range(measurements_per_test):
                frac = (1 - immune_ar[i]) / 2.
                self.grid = np.random.choice([0, 1, 3], size=(self.N, self.M), p=[frac, frac, immune_ar[i]])
                self.update(1, sweeps_per_test, anim=False)
                test_results[j] = self.measure_infected()
            ys[i] = np.average(test_results)
            yerr[i] = np.std(test_results)#  / np.sqrt(measurements_per_test)
        print(yerr)
        if show:
            plt.errorbar(immune_ar, ys, yerr=yerr, capsize=3)
            plt.xlabel("Immune fraction $(f_{Im})$")
            plt.ylabel("Infected fraction")
            plt.title("Immune test")
            plt.show()
        if save:
            np.savetxt('sirs_immune_test.txt', ys)

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
        print("sirs.py M N P1 P2 P3 Pimmune anim:0 or phase:1 or slice:2 or contour:3 or immune:4 or 5:waves")
        sys.exit()
    M = int(argv[0])
    N = int(argv[1])
    P1 = float(argv[2])
    P2 = float(argv[3])
    P3 = float(argv[4])
    P4 = float(argv[5])

    if argv[6] == '0' or argv[6] == 'anim':
        s = Sirs(M, N, p1=P1, p2=P2, p3=P3, immune=P4)
        s.run_animation()
    elif argv[6] == '1' or argv[6] == 'phase':
        s = Sirs(M, N, p1=P1, p2=0.5, p3=P3, immune=0, test=True)
        s.phase_test()
    elif argv[6] == '2' or argv[6] == 'slice':
        s = Sirs(M, N, p1=0.2, p2=0.5, p3=0.5, immune=P4, test=True)
        s.slice_test()
    elif argv[6] == '3' or argv[6] == 'contour':
        s = Sirs(M, N, p1=0.5, p2=0.5, p3=0.5, immune=0, test=True)
        s.contour_test()
    elif argv[6] == '4' or argv[6] == 'immune':
        s = Sirs(M, N, p1=0.5, p2=0.5, p3=0.5, immune=P4, test=True)
        s.immune_test()
    elif argv[6] == '5' or argv[6] == 'waves':
        s = Sirs(M, N, p1=P1, p2=P2, p3=P3, immune=P4, test=True)
        s.waves_test()

main(sys.argv[1:])
