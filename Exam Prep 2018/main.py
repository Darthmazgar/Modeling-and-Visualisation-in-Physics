from fisher_eqn import FisherEqn
import sys

def main(argv):
    if len(argv) != 3:
        print("Start properly.")
        print("python main.py N R dt")
        sys.exit()

    N = int(argv[0])
    R = float(argv[1])
    dt = float(argv[2])

    # dt = 1

    f = FisherEqn(N, R, dt)
    # print(f.phi_grid)
    f.run_animation()


if __name__ == '__main__':
    main(sys.argv[1:])
