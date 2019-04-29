import sys
from Potts import PottsModel

def main(argv):
    if len(argv) != 2:
        print("main.py N T")
        sys.exit()
    N = int(argv[0])
    T = int(argv[0])

    potts = PottsModel(N, T)
    potts.run_animation()






main(sys.argv[1:])
