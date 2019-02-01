from ising import Grid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

def main(argv):
    if len(argv) != 5:
        print("ising_exe.py N M Temp DynamicSystem(0: Glauber, 1: Kawasaki) ('anim':0 / 'test':1, 'run':2)")
        sys.exit()

    N = int(argv[0])
    M = int(argv[1])
    T = float(argv[2])
    ds = int(argv[3])

    if argv[4] == 'anim' or argv[4] == '0':
        grid = Grid(N, M, T, ds=ds, sv_ext='', anim=True, all_up=False)
        # grid.init_kaw_grid()
        grid.animate()
    elif argv[4] == 'test' or argv[4] == '1':
        if ds == 0:
            sv_ext = '_gla'
        elif ds == 1:
            sv_ext = '_kaw'
        else:
            sv_ext = ''
        if ds == 0:
            grid = Grid(N, M, T, ds=ds, sv_ext=sv_ext, anim=False, all_up=True)
        elif ds == 1:
            grid = Grid(N, M, T, ds=ds, sv_ext=sv_ext, anim=False, all_up=False)
            grid.init_kaw_grid()
        grid.temperature_tests()
        grid.susceptibility()
        grid.heat_cap()
        grid.bootstarap_errors()
    elif argv[4] == 'run' or argv[4] == '2':
        grid = Grid(N, M, T, ds=ds, sv_ext='', anim=True, all_up=False)
        sweeps = int(input("Run for how many sweeps: "))
        for i in range(sweeps):
            grid.update_sweep(1)
        grid.imshow_grid()
        print("System magnetisation: %d" % grid.sys_magnetisation())
        print("System energy: %d" % grid.sys_energy())
        plt.show()


main(sys.argv[1:])
