#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
import numpy as np


data11 = np.genfromtxt('./propiedadesVSmasa', comments='#')[:]
radio = data11[:,2]
temperatura = data11[:,1]
luminosidad = data11[:,3]
masa = data11[:,0]
error = 100.0*data11[:,4]


# Plot PRESSURE plot
f, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2, 2, sharex='col')#, sharey='row',

#plt.subplots_adjust(left=0.10, bottom=0.5, right=0.45, top=0.95,wspace=0.40, hspace=0.40)
#ax1 = fig1.add_axes([0.1, 0.14, 0.8, 0.77])
ax1.plot(masa, temperatura ,color="black", lw=1, ls='-',marker='o', markersize=5.1)
ax1.grid()
ax2.plot(masa, radio ,color="black", lw=1, ls='-',marker='o', markersize=5.1)
ax2.grid()
ax3.plot(masa, luminosidad ,color="black", lw=1, ls='-',marker='o', markersize=5.1)
ax3.grid()
ax4.plot(masa, error ,color="black", lw=1, ls='-',marker='o', markersize=5.1)
ax4.grid()



ax3.set_xlabel(r' $ Masa \; \; [10^{33}\; g]$', fontsize=16)
ax4.set_xlabel(r' $ Masa \; \; [10^{33}\; g]$', fontsize=16)

ax1.set_ylabel(r' $ T_c \; \; [10^{7}\;K]$', fontsize=16)
ax2.set_ylabel(r' $ R \; \; [10^{10}\;cm]$', fontsize=16)
ax3.set_ylabel(r' $ L \; \; [10^{33}\;erg/s]$', fontsize=16)
ax4.set_ylabel(r' $ Error\;$%', fontsize=16)


plt.show()



