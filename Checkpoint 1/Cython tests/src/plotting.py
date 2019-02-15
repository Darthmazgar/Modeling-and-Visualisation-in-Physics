import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('magnetisation_kaw.txt')
sus = np.genfromtxt('susceptibility_kaw.txt')
temp = np.genfromtxt('temperature_kaw.txt')
en = np.genfromtxt('energy_kaw.txt')
c = np.genfromtxt('heat_cap_kaw.txt')
e_error = np.genfromtxt('sigma_bs_kaw.txt')


# data = np.genfromtxt('magnetisation_gla.txt')
# sus = np.genfromtxt('susceptibility_gla.txt')
# temp = np.genfromtxt('temperature_gla.txt')
# en = np.genfromtxt('energy_gla.txt')
# c = np.genfromtxt('heat_cap_gla.txt')
# e_error = np.genfromtxt('sigma_bs_gla.txt')

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
plt.errorbar(temp, c, yerr=e_error, capsize=5)
# plt.plot(temp, c)
plt.title("Heat Capacity vs Temperature")
plt.ylabel("Heat Capacity")
plt.xlabel("Temperature")
plt.show()
