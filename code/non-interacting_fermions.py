# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 11:45:18 2021

@author: roshn
"""
import numpy as np 
import matplotlib.pyplot as plt 
from numpy.linalg import eig
from mpl_toolkits.mplot3d import Axes3D

# 2 identical non-interacting fermions (spin-1/2 fermions) in an infinite potential well 
#modelling the eigenfunctions and energy levels 

#note that identical means - they cannot be distinguished, i.e. their spins are parallel 

#Solving the TISE for 1 particle in an infinite potential well using finite difference method 

#initial constants 
hbar = 1.05*(10**-34)
m = 9.1*(10**-31)      #electron mass in atomic units 


#well width 
W = (10**-10)

#step size  
n = 1000

#x-axis points and y-axis points 
x = np.linspace(0,W,n)
y = np.linspace(0,W,n)

#print(x)
#number of grid points 
h = x[1] - x[0]



#Kinetic energy matrix 
K = np.zeros((n,n))
s,t = np.shape(K)
for i in range(s):
    for j in range(t):
            if i==j: 
                K[i,j] = -2 
for i in range(s):
    for j in range(t-1):
            if j==(i+1):
                K[i,j] = 1 
for i in range (1,s):
    for j in range(t):
            if j==(i-1):
                K[i,j] = 1 

#print(K)

                
#Potential Energy function 
def P(x):
    return 0 
            
#Potential Energy Matrix     
V = np.zeros((n,n))
for i in range(s):
    for j in range(t):
        if i==j: 
            V[i,j] = P(x)
        else: continue 

#print(V)
        
#Hamiltonian 
        
H = (-hbar**2/(2*m*(h**2)))*K + V 

#Solving eigenvalue problem 
eigval, eigvec = eig(H)


#sorting eigenvalues in asceding order and arranging eigenvectors too 
idx = np.argsort(eigval)
eigval = eigval[idx]
eigvec = eigvec[:,idx]

#Plotting 
level1 = int(input('Enter energy level of electron 1:'))
level2 = int(input('Enter energy level of electron 2:'))


#mesh points for plotting in 3-D 
xx, yy = np.meshgrid(x,y)  

if level1 != level2: 
    
    #Wavefunction matrix  
    Z = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Z[i,j] = (eigvec[i:i+1,level1-1]*eigvec[j:j+1,level2-1] - eigvec[j:j+1,level1-1]*eigvec[i:i+1,level2-1])
        


    #Probability Density matrix  
    Z2 = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Z2[i,j] = np.abs(Z[i,j])**2
    
    #Normalization constant 
    A = np.trapz(np.trapz(Z2,x),y)
    
    Z = (1/np.sqrt(A))*Z
    Z2 = (1/A)*Z2
    
    fig= plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    p = ax.plot_surface(xx,yy,Z, cmap='hot')
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title(' Normalized Wavefunction')
        
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    p = ax.plot_surface(xx,yy,Z2, cmap='hot')
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title('Normalized Probability Density')
    
    E_total = (eigval[level1-1] + eigval[level2-1])*(6.242*(10**18))
    print('Total Energy (in eV):',E_total)

else: 
    print("No state possible. Identical fermions cannot occupy the same energy level!")
    Z = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Z[i,j] = 0
    fig= plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    p = ax.plot_surface(xx,yy,Z, cmap='hot')
    ax.set_xlabel("$x_1$",fontsize = '20')
    ax.set_ylabel("$x_2$",fontsize = '20')
    ax.set_title('Normalized Wavefunction', fontsize = '20')
    
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    p = ax.plot_surface(xx,yy,Z, cmap='hot')
    ax.set_xlabel("$x_1$",fontsize = '20')
    ax.set_ylabel("$x_2$",fontsize = '20')
    ax.set_title('Normalized Probability Density', fontsize = '20')


    
                                                #ANALYTICAL
def psix(x,l,a):
    return np.sin((l*np.pi*x)/a)
def psiy(y,l,a):
    return np.sin((l*np.pi*y)/a)

def PSI(x,y,a,l1,l2):
    return (np.sqrt(2)/a)*(psix(x, l1, a)*psiy(y, l2, a) - psix(x, l2, a)*psiy(y, l1, a))

#fig = plt.figure()
#ax = fig.add_subplot(1, 2, 1, projection='3d')
#p = ax.plot_surface(xx,yy,PSI(xx,yy,W,level1,level2), cmap='cool')
    

                                                    #***
    
#Perturbation: considering the potential due to the interaction between 2 electrons 
    
e = 1.6*(10**-19)
eo = 8.85*(10**-12)

                                
                            #FIRST ORDER CORRECTION OF ENERGY
        
P = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i != j:
            P[i,j] = Z2[i,j]*((e**2)/(4*np.pi*eo*np.abs(i-j)))
        else: 
            P[i,j] = Z2[i,j]*0
            
#print(P)

E_correction = (np.trapz(np.trapz(P,x),y))*(6.242*(10**18))

E_corrected = (E_correction + E_total)
print('Corrected Energy (in eV):', E_corrected)
print("Correction in energy (in eV):",E_correction)




