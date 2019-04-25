import sys
from diffusion import *

def main(argv):
    if len(argv) != 3:
        print("main.py N kapa sigma")
        sys.exit()
    N = int(argv[0])
    kapa = float(argv[1])
    sigma = float(argv[2])


    dif = Diffusion(N, kapa=kapa, sigma=sigma)
    dif.run_animation()



    # dif.animation = False
    # dif.measure_mean_psi()
    # dif.psi_with_r()




main(sys.argv[1:])
