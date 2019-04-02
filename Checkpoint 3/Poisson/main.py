import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys
from poisson import *

def main(args):
    if len(args) != 1:
        print("\nTo few arguments.\n")
        print("python main.py N")
        sys.exit()
    N = int(args[0])
    # dt = float(args[1])
    # dx = float(args[2])
    # phi_0 = float(args[3])
    # test = args[4]

    p = Poisson(N, 0.1)
    # print(p.grid)
    # p.add_point_charge(int(N/2), int(N/2), int(N/2))
    p.add_line_charge(int(N/2), int(N/2))
    # p.run_animation()
    # p.update(1)
    # print(p.phi_grid)
    print(p.rho_grid)


    # p.run_animation()



if __name__ == '__main__':
    main(sys.argv[1:])
