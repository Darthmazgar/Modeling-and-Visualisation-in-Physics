import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys
from poisson import *

def main(args):
    if len(args) != 4:
        print("\nTo few arguments.\n")
        print("python main.py N (anim or plot or test) method(jacobi, gauss_seidel) distribution(line or point)")
        sys.exit()
    N = int(args[0])
    method = args[2]
    p = Poisson(N, 0.1, method=method)
    if args[3] == 'point':
        p.add_point_charge(int(N/2), int(N/2), int(N/2), val=1.)
    if args[3] == 'line':
        p.add_line_charge(int(N/2), int(N/2), val=1.)
    if args[3] == 'dipole':
        p.add_point_charge(int(N/3), int(N/3), int(N/3), val=1.)
        p.add_point_charge(int(2*N/3), int(2*N/3), int(2*N/3), val=-1.)
    if args[3] == 'plines':
        p.add_line_charge(int(N/3), int(N/3), val=1.)
        p.add_line_charge(int(2*N/3), int(2*N/3), val=1.)

    p.sweeps = 1000

    if args[1] == 'anim':
        p.sweeps = 1
        p.run_animation()
    if args[1] == 'test':
        p.contour_test()
        sys.exit()
    else:
        p.update(1)
        p.plot_E_field(save=True)
        p.plot_B_field(save=True)

if __name__ == '__main__':
    main(sys.argv[1:])
