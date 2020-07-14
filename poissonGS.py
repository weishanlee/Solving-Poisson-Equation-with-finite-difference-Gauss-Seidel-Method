# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:21:15 2020
@author: Wei-shan
Solving Poisson Equation 
    del^2 V = -rho
in the electric potential problem 
with Finite difference of Gauss-Seidel Method and overrelaxation.
Reference: Mark Newman, Computational Physics, CH9.
"""
from pylab import imshow,gray,show
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd

# Constants
M = 100         # Grid squares on a side
V = 30.0         # Voltage at top wall
rho = 1.0       # charge density
target = 1e-6   # Target accuracy
omega = 0.9
delta = 1.0
xMin = 0.0
xMax = 100.0
yMin = 0.0
yMax = 100.0

# Create arrays to hold potential values
phi = np.zeros([M+1,M+1],float) # phi[y,x]

# Main loop
while delta>target:
    delta = 0.0
    # Calculate new values of the potential
    for i in range(M+1):
        for j in range(M+1):
            ## Setting up boundary conditions (Simulation 1)
            if i==0:
                phi[i,j] = V
            elif j==0:
                phi[i,j] = V
            elif i==M:
                phi[i,j] = -V
            elif j==M:
                #or j==M: # something is wrong for setting phi[i,j]=V at i==M or j==M
                phi[i,j] = -V
            ## End of Setting up boundary conditions (Simulation 1No)
            ## Setting up other boundary values inside the digram (Simulation 3)
            #if i==0 or i==M or j==0 or j==M:
            #    phi[i,j] = 0.0 
            #elif (i>=20 and i<=80) and j==20:
            #    phi[i,j]=V
            #elif (i>=20 and i<=80) and j==80:
            #    phi[i,j]=-V
            ## End of Setting up other boundary values inside the digram (Simulation 3)
            ## Charge densit
            ## Setting up rho value (Simulation 2)
            elif ( ( i>=20 and i<=40 ) and ( j>=60 and j<=80 ) ):
                temp = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])*(1+omega)/4 - omega * phi[i,j] + 1/4*rho
                if ( abs( phi[i,j] - temp ) > delta ): delta = abs( phi[i,j] - temp )
                phi[i,j] = temp
            #elif ( ( i>=60 and i<=80 ) and  ( j>=20 and j<=40 ) ):
            #    temp = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])*(1+omega)/4 - omega * phi[i,j] - 1/4*rho
            #    if ( abs( phi[i,j] - temp ) > delta ): delta = abs( phi[i,j] - temp )
            #    phi[i,j] = temp
            ## End of Setting up rho value (Simulation 2)
            else:
                temp = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])*(1+omega)/4 - omega * phi[i,j]
                if ( abs( phi[i,j] - temp ) > delta ): delta = abs( phi[i,j] - temp )
                phi[i,j] = temp
# End of main loop

# Make a plot
ax = plt.gca()
imshow(phi,extent=[xMin,xMax,yMin,yMax])#,origin="lower")

plt.minorticks_on()
minorLocatorX = AutoMinorLocator(4) # number of minor intervals per major 
                                    # inteval
minorLocatorY = AutoMinorLocator(4)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
ax.set_xticklabels(ax.get_xticks(),family='monospace',fontsize=10)
ax.set_yticklabels(ax.get_yticks(),family='monospace',fontsize=10)
gray()
show()

# Save location vs voltage into csv file.
phiData = pd.DataFrame(columns = ['X','Y','Voltage'])
phiData_file = open(r'E:\github\Solving-Poisson-Equation-with-finite-difference-Gauss-Seidel-Method\poissonGS.csv','w',newline='') 

X = []
Y = []
Voltage = []

for i in range(M+1):
    for j in range(M+1):
        X += [j * (xMax-xMin)/M]
        Y += [i * (yMax-yMin)/M]
        Voltage += [ phi[i,j] ]

phiData['X'] = X 
phiData['Y'] = Y 
phiData['Voltage'] = Voltage

phiData.to_csv(phiData_file, sep=',', encoding='utf-8', index=False) 
phiData_file.close()