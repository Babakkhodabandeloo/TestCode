import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,2,0.01)
f = 1 # Hz
y = 2*np.sin(2*np.pi*f*x) + 0.5*np.cos(2*np.pi*3*f*x)

plt.plot(x,y)
plt.show()
