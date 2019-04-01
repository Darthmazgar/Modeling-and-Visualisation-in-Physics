import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('magnetisation.txt')
sus = np.genfromtxt('susceptibility.txt')
temp = np.genfromtxt('temperature.txt')
en = np.genfromtxt('energy.txt')
c = np.genfromtxt('heat_cap.txt')
magnetisation = [np.average(data[x]) for x in range(len(data))]
energy = [np.average(en[x]) for x in range(len(en))]

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
plt.plot(temp, energy)
plt.title("Energy vs Temperature")
plt.ylabel("Energy")
plt.xlabel("Temperature")
plt.show()
plt.plot(temp, c)
plt.title("Heat Capacity vs Temperature")
plt.ylabel("Heat Capacity")
plt.xlabel("Temperature")
plt.show()
