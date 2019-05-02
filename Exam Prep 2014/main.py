import sys
from superfluid import SuperfluidTransiton

def main(argv):
    if len(argv) != 4:
        print("main.py L T c anim/test")
        sys.exit()
    L = int(argv[0])
    T = float(argv[1])
    c = float(argv[2])
    if c > 1 or c <= 0:
        c = float(input("c is the fraction of sites which are He4 so muct be between 0 and 1:\n"))
    if L <= 0:
        L = int(input("L must be a positive integer:\n"))

    i = SuperfluidTransiton(L, T, c)
    if argv[3] == 'anim':
        i.run_animation()
    elif argv[3] == 'test':
        i.run_tests()


if __name__ == '__main__':
    main(sys.argv[1:])
