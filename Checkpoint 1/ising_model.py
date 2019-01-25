import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class DynamicSystem:
    def __init__(self, T):
        self.T = T
        self.J = 1
        self.kb = 1  # 1.38*10**-23

    def set_T(self, T):
        self.T = T
        return self.T

    def get_T(self):
        return self.T

    def get_Kb(self):
        return self.kb

    def calc_energy(self, points):
        sample = points[0]
        energy = 0
        for i in range(1, len(points)):
            energy += -self.J * (sample * points[i])
        return energy

    def P(self, delta_E):
        if self.T == 0:
            return 0
        return np.exp(- (delta_E / (self.kb * self.T)))

    def glauber_dynamics(self, points, coords):
        """
        Takes an array of one point and calculates the energy change and if it
        is favourable.
        """
        init_E = self.calc_energy(points[0])
        points[0][0] *= -1  # Switch state of choice point.
        final_E = self.calc_energy(points[0])
        delta_E = final_E - init_E
        if delta_E <= 0:
            return True  # If true then flip
        else:
            rand = np.random.random()
            prob = self.P(delta_E)
            if rand <= prob:
                return True
            else:
                return False

    def kawasaki_dynamics(self, points, coords):
        # TODO Consider double counting for nearest neighbours
        point1 = points[0]
        point2 = points[1]
        if point1[0] == point2[0]:
            return False
            # pass  # Return whatever happens if there is no change.
        else:
            init_E = self.calc_energy(point1) + self.calc_energy(point2)
            point1[0] *= -1
            point2[0] *= -1
            final_E = self.calc_energy(point1) + self.calc_energy(point2)
            delta_E = final_E - init_E
            if delta_E <= 0:  # The following is assuming that there are the same steps as for grober
                return True  # If true then flip
            else:
                rand = np.random.random()
                prob = self.P(delta_E)
                if rand <= prob:
                    return True
                else:
                    return False


class Lattice:
    def __init__(self, N, sys, points_needed, ds, sps=2500, all_up=True):
        self.N = N
        self.sys = sys
        self.pneeded = points_needed  # Number of random points needed
        self.ds = ds
        self.steps_per_sweep = sps
        self.grid = np.ones((N, N), dtype=int)
        if not all_up:
            self.init_grid()
        self.fig = plt.figure()

    def print_grid(self):
        print(self.grid)

    def imshow_grid(self):
        plt.title("Ising Model")
        plt.imshow(self.grid, cmap='Blues', interpolation='sinc', vmin=-1, vmax=1)
        plt.colorbar()

    def init_grid(self):
        spin = np.array([-1, 1])
        for j in range(len(self.grid)):
            for i in range(len(self.grid[0])):
                self.grid[i][j] = np.random.choice(spin)
        return self.grid

    def select_rand_coords(self):
        sel = np.random.randint(self.N, size=2)
        return sel

    def around_selection(self, sel):
        """
        Creates array of values of nearest neighbours around a selection.
        :param sel: (np.array(2)), Selected point in self.grid.
        :return: (np.array(5)), [Selected val, one down, one up, one right, one left]
        """
        ar = np.zeros(5)
        ar[0] = self.grid[sel[0]][sel[1]]
        if sel[0] == self.N - 1:
            sel[0] = -1
        if sel[1] == self.N - 1:
            sel[1] = -1
        ar[1] = self.grid[sel[0] + 1][sel[1]]
        ar[2] = self.grid[sel[0] - 1][sel[1]]
        ar[3] = self.grid[sel[0]][sel[1] + 1]
        ar[4] = self.grid[sel[0]][sel[1] - 1]
        return ar

    def test_selection(self, sel):
        points_around = self.around_selection(sel)
        if self.sys.favourable(points_around):
            if self.grid[sel[0]][sel[1]] == 0:
                self.grid[sel[0]][sel[1]] = 1
            else:
                self.grid[sel[0]][sel[1]] = 0

    def save_grid(self, outfile):
        np.savetxt(outfile, self.grid)

    def run_sim(self, sweeps, save=False):
        for i in range(self.steps_per_sweep * sweeps):
            selections = []
            points = []
            for j in range(self.pneeded):
                selections.append(self.select_rand_coords())
                points.append(self.around_selection(selections[j]))

            if self.sys(points, selections):
                for k in range(len(selections)):
                    self.grid[selections[k][0]][selections[k][1]] *= -1
            if save:
                self.save_grid("grid.txt")

    def animate_update(self, i, save=False):
        for i in range(self.steps_per_sweep):
            selections = []
            points = []
            for j in range(self.pneeded):
                selections.append(self.select_rand_coords())
                points.append(self.around_selection(selections[j]))

            if self.sys(points, selections):
                for k in range(len(selections)):
                    self.grid[selections[k][0]][selections[k][1]] *= -1
            if save:
                self.save_grid("grid.txt")

        self.fig.clear()
        self.imshow_grid()

    def animate(self):
        anim = FuncAnimation(self.fig, self.animate_update)
        plt.show()

    def sys_magnetisation(self):
        """
        Magnetisation per data point in the simulation
        :return: <|m|> = |m_i/n|
        """
        M = np.sum(self.grid)
        mag = np.abs(M )
        return mag

    def sys_energy(self):
        E = 0
        for j in range(len(self.grid)):
            for i in range(len(self.grid[j])):
                point = self.around_selection([i, j])
                E += self.ds.calc_energy(point)
        E /= 2  # Prevent double counting.
        return E

    def temperature_tests(self, t_min=1, t_max=3, data_points=20, sweeps=100, tests=1000, eng=True, mag=True, save=True):
        temperature = np.linspace(t_min, t_max, data_points)
        magnetisation = np.zeros((data_points, tests))
        energy = np.zeros((data_points, tests))

        for i in range(data_points):
            sys.stdout.write("Simulation progress: %.1f%%\r" % (100 * i / data_points))
            sys.stdout.flush()

            self.ds.set_T(temperature[i])  # Set the temperature of the system.
            self.run_sim(sweeps)
            for j in range(tests):
                self.run_sim(10)
                if mag:
                    magnetisation[i][j] = self.sys_magnetisation()
                if eng:
                    energy[i][j] = self.sys_energy()

        if save:
            np.savetxt('temperature_kaw.txt', temperature)
            if mag:
                np.savetxt('magnetisation_kaw.txt', magnetisation)
            if eng:
                np.savetxt('energy_kaw.txt', energy)

    def susceptibility(self, save=True):
        data = np.genfromtxt('magnetisation_kaw.txt')
        temp = np.genfromtxt('temperature_kaw.txt')
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        chi = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N**2 * self.ds.kb * temp[i])
            chi[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt('susceptibility_kaw.txt', chi)

    def heat_cap(self, save=True):
        data = np.genfromtxt('energy_kaw.txt')
        temp = np.genfromtxt('temperature_kaw.txt')
        magnetisation = [np.average(data[x]) for x in range(len(data))]
        C = np.zeros(len(temp))
        for i in range(len(temp)):
            norm_fact = 1 / (self.N**2 * self.ds.kb * temp[i]**2)
            C[i] = norm_fact * (np.average(np.square(data[i])) - np.square(np.average(data[i])))
        if save:
            np.savetxt('heat_cap_kaw.txt', C)


def main():
    # T = float(input("Enter the temperature of the system: "))
    ds = DynamicSystem(2.3)

    # dynamic_sys = [(ds.glauber_dynamics, 1), (ds.kawasaki_dynamics, 2)]
    # sys = int(input("Choose the system dynamics: 0, %s: 1, %s: " % (dynamic_sys[0], dynamic_sys[1])))
    # print(dynamic_sys[sys][0])

    lattice = Lattice(50, ds.glauber_dynamics, 1, ds=ds, all_up=True)
    # lattice = Lattice(50, ds.kawasaki_dynamics, 2, ds=ds, all_up=False)

    # lattice.run_sim(100)  # Run for a certain number of sweeps.
    # lattice.imshow_grid()  # Display grid after n sweeps.
    # lattice.sys_magnetisation()
    # plt.show()
    #
    lattice.animate()  # Animate live
    #
    # lattice.temperature_tests()  # Run Tests
    # lattice.susceptibility()
    # lattice.heat_cap()
    #
    # lattice.sys_magnetisation()
    # lattice.temperature_test()

main()
