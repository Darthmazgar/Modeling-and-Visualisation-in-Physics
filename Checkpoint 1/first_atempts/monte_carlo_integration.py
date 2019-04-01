import numpy as np
import matplotlib.pyplot as plt
import itertools

class McInt:
    def __init__(self, N, fn):
        self.N = N
        self.fn = fn
        self.ar = False

    def integrate(self, plot=True):
        count = 0
        col = (1, 0, 0)
        ar = np.zeros((self.N, 3))
        for i in range(self.N):
            rand_x = np.random.random()*2 - 1
            rand_y = np.random.random()*2 - 1
            ar[i] = np.array([rand_x, rand_y, False])
            if self.fn(rand_x, rand_y):
                count += 1
                col = (0, 1, 0)
                ar[i][2] = True
        print((count / self.N))
        self.ar = ar
        return count / self.N

    def plot(self):
        BLUE = (0, 0, 1)
        RED = (1, 0, 0)
        xs = [self.ar[i][0] for i in range(len(self.ar))]
        ys = [self.ar[i][1] for i in range(len(self.ar))]
        col = [RED if self.ar[i][2] == 0 else BLUE for i in range(len(self.ar))]
        print(col)
        plt.plot(xs, ys, 'x', color=col)
        plt.show()


class Function:
    @staticmethod
    def x_squared(x, y):
        if y <= x**2:
            return True
        else:
            return False

    @staticmethod
    def x(x, y):
        if y <= x:
            return True
        else:
            return False

    @staticmethod
    def circle(x, y):
        if np.abs(y) <= 1 - (x**2 + y**2) and np.abs(x) <= 1 - (x**2 + y**2):
            return True
        else:
            return False


def main():
    mc = McInt(1000, Function.circle)
    mc.integrate()
    # mc.plot()


main()
