"""
2019 Modeling and Visualisation in Physics Exam
B075090
"""

import sys
from autocat import Autocat

def main(argv):
    if len(argv) != 5:
        print("python3 main.py N R F dt anim/testc/testd")
        sys.exit()
    N = int(argv[0])
    R = float(argv[1])
    F = float(argv[2])
    dt = float(argv[3])

    i = Autocat(N, R, F, dt=dt)
    if argv[4] == 'anim':
        i.run_animation()
    elif argv[4] == 'testc':
        i.run_c_tests()
    elif argv[4] == 'testd':
        i.run_d_tests()


if __name__ == '__main__':
    main(sys.argv[1:])
