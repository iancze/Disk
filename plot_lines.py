import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
img = np.load("img45.npy")
plt.plot(np.arange(0,50,1),img[:,24])
plt.plot(np.arange(0,50,1),img[:,18])
plt.show()
