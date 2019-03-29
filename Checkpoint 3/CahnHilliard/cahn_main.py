import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import sys
from cahn_hilliard import *


def main(args):
    if len(args) != 5:
        print("\nTo few arguments.\n")
        print("python cahn_hilliard.py N dt dx phi_0 0:anim or 1:test or 2:run")
        sys.exit()
    N = int(args[0])
    dt = float(args[1])
    dx = float(args[2])
    phi_0 = float(args[3])
    test = args[4]

    if args[4] == '0' or args[4] == 'anim':
        ch = CahnHil(N, dt=dt, dx=dx, phi_0=phi_0, anim='True')
        if len(args) == 6:
            for i in range(int(args[5])):
                ch.add_drop(np.random.uniform()*15, np.random.randint(N),
                        np.random.randint(N))
        ch.steps_per_update = 100
        ch.run_animation()
    elif args[4] == '1' or args[4] == 'test':
        ch = CahnHil(N, dt=dt, dx=dx, f_test=True,phi_0=phi_0)
        ch.steps_per_update = 200
        ch.run_f_test(100)
    elif args[4] == '2' or args[4] == 'run':
        sweeps = int(input("Run for how many sweeps?\n"))
        ch = CahnHil(N, dt=dt, dx=dx, f_test=True,phi_0=phi_0, anim='True')
        ch.steps_per_update = sweeps
        ch.update(1, anim=args[4])
        # cd.display(sweeps)
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
