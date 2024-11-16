import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,2,0.01)
f1 = 1 # Hz
f2 = 4 # Hz
y = 2*np.sin(2*np.pi*f1*x) + 0.3*np.cos(2*np.pi*f2*x)

plt.plot(x,y)
plt.show()
