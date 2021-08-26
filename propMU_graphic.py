#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
import numpy as np


data11 = np.genfromtxt('./propiedadesVSmetal', comments='#')[:]
radio1 = data11[0:6,5]
temperatura1 = data11[0:6,3]
luminosidad1 = data11[0:6,4]
metal1 = 100*data11[0:6,2]
error1 = data11[0:6,6]

radio2 = data11[7:13,5]
temperatura2 = data11[7:13,3]
luminosidad2 = data11[7:13,4]
metal2 = 100*data11[7:13,2]
error2 = data11[7:13,6]

mu1 = 1.0/(2.0*data11[0:6,0] + 3.0/4.0*data11[0:6,1] + 0.5*data11[0:6,2])
mu2 = 1.0/(2.0*data11[7:13,0] + 3.0/4.0*data11[7:13,1] + 0.5*data11[7:13,2])

# Plot PRESSURE plot
f, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2, 2, sharex='col')#, sharey='row',

#plt.subplots_adjust(left=0.10, bottom=0.5, right=0.45, top=0.95,wspace=0.40, hspace=0.40)
#ax1 = fig1.add_axes([0.1, 0.14, 0.8, 0.77])
ax1.plot(mu1, temperatura1 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando Y')
ax1.grid()
ax2.plot(mu1, radio1 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando Y')
ax2.grid()
ax3.plot(mu1, luminosidad1 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando Y')
ax3.grid()
ax4.plot(mu1, error1 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando Y')
ax4.grid()

ax1.plot(mu2, temperatura2 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando X')
ax1.grid()
ax2.plot(mu2, radio2 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando X')
ax2.grid()
ax3.plot(mu2, luminosidad2 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando X')
ax3.grid()
ax4.plot(mu2, error2 ,color="purple", lw=0, ls='-',marker='o', markersize=5.1, label ='Variando X')
ax4.grid()
"""
ax1.legend(loc=0,prop={'size':9})
ax2.legend(loc=0,prop={'size':9})
ax3.legend(loc=0,prop={'size':9})
ax4.legend(loc=0,prop={'size':9})
"""


ax3.set_xlabel(r' $ \mu $', fontsize=16)
ax4.set_xlabel(r' $ \mu $', fontsize=16)

ax1.set_ylabel(r' $ T_c \; \; [10^{7}\;K]$', fontsize=16)
ax2.set_ylabel(r' $ R \; \; [10^{10}\;cm]$', fontsize=16)
ax3.set_ylabel(r' $ L \; \; [10^{33}\;erg/s]$', fontsize=16)
ax4.set_ylabel(r' $ Error\;$ %', fontsize=16)


plt.show()



