from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import stats
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pylab import *
from uncertainties import ufloat
from mpl_toolkits.mplot3d import Axes3D

plt.clf()

fig = plt.figure(figsize=(8,6))

ax = fig.add_subplot(111,projection='3d')

x,y,z,N=np.loadtxt("Nachbarnt0.txt", unpack=True)

colors = cm.hsv(N/max(N))

colmap = cm.ScalarMappable(cmap=cm.hsv)
colmap.set_array(N)

yg = ax.scatter(x, y, z, c=colors, marker='o')
cb = fig.colorbar(colmap)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

