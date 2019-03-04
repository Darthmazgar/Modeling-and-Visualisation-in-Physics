from ising import Grid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import timeit
# from libcpp cimport bool


# # py = timeit.timeit('cy_ising_model.main()', setup='import cy_ising_model', number=1)
# cy = timeit.timeit('ising.Grid.main()', setup='import ising', number=1)
#
# print(cy, py)
# print("Cython is {}x faster".fotmat(py/cy))


def main():
    grid = Grid(50, 50, 2, anim=True, all_up=False)
    # grid.print_grid()
    # grid.imshow_grid()
    # plt.show()
    # for i in range(100):
    #     grid.update_sweep(1)
    # # grid.print_grid()
    # grid.imshow_grid()
    # plt.show()
    # grid.animate()
    # grid.temperature_tests()
    # grid.susceptibility()
    # grid.heat_cap()
    print(grid.nn_check(0, 0, 0, 1))

main()
