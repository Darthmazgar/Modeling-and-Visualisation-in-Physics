import sys
from other import Other ##################

def main(argv):
    if len(argv) != 2:
        print("main.py N anim/test")
        sys.exit()
    N = int(argv[0])

    i = Other(N) ##############
    if argv[1] == 'anim':
        i.run_animation()
    elif argv[1] == 'test':
        i.run_tests()


if __name__ == '__main__':
    main(sys.argv[1:])
