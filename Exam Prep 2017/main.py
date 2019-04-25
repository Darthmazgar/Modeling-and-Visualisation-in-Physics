import sys
from ca import *

def main(argv):
    if len(argv) != 3:
        print("main.py N P1 P2 ")
        sys.exit()
    N = int(argv[0])
    p1 = float(argv[1])
    p2 = float(argv[2])

    ca = CellularAutomata(N, p1, p2)
    ca.run_animation()
    # ca.animation = False
    # ca.vary_p2()




main(sys.argv[1:])
