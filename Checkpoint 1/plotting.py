import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('Magnetisation_kaw.txt')
sus = np.genfromtxt('susceptibility_kaw.txt')
temp = np.genfromtxt('temperature_kaw.txt')
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
