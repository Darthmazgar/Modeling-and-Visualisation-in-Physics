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

    p = Poisson(N, 0.1)
    # print(p.grid)
    p.add_point_charge(int(N/2), int(N/2), int(N/2))
    # p.run_animation()
    # p.update(1, 10)
    p.update(1)
    p.plot_E_field()

if __name__ == '__main__':
    main(sys.argv[1:])
