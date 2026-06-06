import matplotlib.pyplot as plt
from matplotlib import gridspec

fig = plt.figure()

gs = fig.add_gridspec(2, 4)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 1:3])
ax3 = fig.add_subplot(gs[0, 2:4])
#ax3 = fig.add_subplot(gs[:, 3])

plt.show()