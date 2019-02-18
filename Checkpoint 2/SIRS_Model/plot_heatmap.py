import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

grid = np.genfromtxt('sirs_heatmap.txt')
print(np.max(grid))
plt.imshow(grid, interpolation='bilinear',
           cmap='Blues', vmin=0, vmax=0.5)
plt.xlabel("P3")
plt.ylabel("P1")
plt.title("Sirs test")
# plt.colorbar()
plt.show()
