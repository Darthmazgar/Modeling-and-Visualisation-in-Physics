import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('magnetisation.txt')
sus = np.genfromtxt('susceptibility.txt')
temp = np.genfromtxt('temperature.txt')
c = np.genfromtxt('heat_cap.txt')
magnetisation = [np.average(data[x]) for x in range(len(data))]
plt.plot(temp, magnetisation)
plt.title("Magnetisation vs Temperature")
plt.ylabel("Magnetisation")
plt.xlabel("Temperature")
plt.show()
plt.plot(temp, sus)
plt.title("Susceptibility vs Temperature")
plt.ylabel("Susceptibility")
plt.xlabel("Temperature")
plt.show()
plt.plot(temp, c)
plt.title("Heat Capacity vs Temperature")
plt.ylabel("Heat Capacity")
plt.xlabel("Temperature")
plt.show()
