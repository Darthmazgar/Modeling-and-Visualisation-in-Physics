from ising import *
import sys


def main(argv):
    if len(argv) != 1:
        print("python main.py N")
    N = int(argv[0])

    i = Ising(N=N)
    i.run_animation()


if __name__ == '__main__':
    main(sys.argv[1:])
