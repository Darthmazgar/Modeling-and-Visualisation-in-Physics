import numpy as np
import matplotlib.pyplot as plt

def get_g(vals):
    return np.average(np.gradient((vals)))

data = np.genfromtxt('output')
xs = data[:, 0]
ys = data[:, 1]
x_grad = get_g(xs)
y_grad = get_g(ys)
m_g = np.sqrt(x_grad**2 + y_grad**2)
print("Gradient: %.4f" % m_g)
width = int(25)
plt.plot(xs, ys)
plt.text(width- 7, 0, "Gradient: %.4f" % m_g)
plt.xlabel("Mean x displacment.")
plt.ylabel("Mean y displacment.")
plt.title("Glider position changing with time.")
plt.show()
