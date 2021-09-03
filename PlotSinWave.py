import matplotlib.pyplot as plt
import numpy as np

x=np.arange(0,1,0.001)
y=np.sin(2*np.pi*4*x)
z=np.cos(2*np.pi*4*x)
plt.plot(x,y)
plt.plot(x,z)