# Angell Plot

import numpy as np
import matplotlib.pyplot as plt


temp_glass = 200
temp_vogel = temp_glass - 50

temp_A = 1727

preF = 34.5

eta_inf = 10**(-3) # unit Pa s

temp_min, temp_max = temp_vogel, 300
temp_max2 = 10000

temp = np.linspace(temp_min, temp_max, 50)
temp2 = np.linspace(temp_min, temp_max2, 1000)
scaled_temp = [temp_glass/T for T in temp]
scaled_temp2 = [temp_glass/T for T in temp2]
Arrhenius = [eta_inf*np.exp(preF*temp_glass/T) for T in temp]
VFT = [eta_inf*np.exp(temp_A/(T - temp_vogel)) for T in temp]
log_Arrhenius = [np.log10(arr) for arr in Arrhenius]
log_VFT = [np.log10(vft) for vft in VFT]
Arrhenius2 = [eta_inf*np.exp(preF*temp_glass/T) for T in temp2]
VFT2 = [eta_inf*np.exp(temp_A/(T - temp_vogel)) for T in temp2]
log_Arrhenius2 = [np.log10(arr) for arr in Arrhenius2]
log_VFT2 = [np.log10(vft) for vft in VFT2]
func_min = np.min(log_Arrhenius)
func_max = np.max(log_Arrhenius)
ylim = [0, func_max]

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title('Arrhnius vs Vogel-Fulcher-Tammann')
ax.set_xlim(0, 1)
ax.set_ylim(-3, 13)
ax.set_xlabel(r'$T_{{g}}$ / $T$')
ax.set_ylabel(r'log($\eta$($T$) /Pa s)')
ax.scatter(scaled_temp, log_Arrhenius, label='Arrhenius')
ax.scatter(scaled_temp, log_VFT, label='Vogel-Fulcher-Tammann')
ax.plot(scaled_temp2, log_Arrhenius2, ls='--', label='strong liquid')
ax.plot(scaled_temp2, log_VFT2, ls='--', label='fragile liquid')

ax.legend(loc='upper left') 
ax.grid()
ax.set_axisbelow(True)

fig.savefig('./png/Angell.png', dpi=300)

plt.show()