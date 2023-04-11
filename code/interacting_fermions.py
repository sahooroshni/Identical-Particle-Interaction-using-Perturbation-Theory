# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 22:57:49 2021

@author: roshn
"""
import numpy as np 
import matplotlib.pyplot as plt 
from numpy.linalg import eig
from mpl_toolkits.mplot3d import Axes3D


#initial constants in SI units 
hbar = 1.05*(10**-34)
m = 9.1*(10**-31)          #electron mass in atomic units 

#well width 
W = 10*(10**-9)             #10 nanometer 

#total number of grid points   
n = 1000

#x-axis points and y-axis points 
x = np.linspace(0,W,n)
y = np.linspace(0,W,n)

#print(x)

#step size 
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
        
#Hamiltonian Matrix 
        
H = (-hbar**2/(2*m*(h**2)))*K + V 

#Solving eigenvalue problem 
eigval, eigvec = eig(H)


#sorting eigenvalues in asceding order and arranging eigenvectors too 
idx = np.argsort(eigval)
eigval = eigval[idx]
eigvec = eigvec[:,idx]

#function that gives the wavefunction and probability density of any energy level 

def Psi(l1,l2):   
    
    #Wavefunction matrix  
    Z = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Z[i,j] = (eigvec[i:i+1,l1-1]*eigvec[j:j+1,l2-1] - eigvec[j:j+1,l1-1]*eigvec[i:i+1,l2-1])
        


    #Probability Density matrix  
    Z2 = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Z2[i,j] = np.abs(Z[i,j])**2
    
  
    
    E_total = (eigval[l1-1] + eigval[l2-1])
    
    return{'psi':Z,'psi2':Z2, 'energy':E_total}


# function that returns the scalar term in the formula for perturbed wavefunctions 
    
def K(l1,l2): 
    
    #well width 
    W = 10*(10**-9)

    #step size  
    n = 1000
    
    #x-axis points and y-axis points 
    x = np.linspace(0,W,n)
    y = np.linspace(0,W,n)
        
    #print(x)
    
    psi1 = Psi(1,3)['psi']        #enter l1 and l2 depending on the energy level at which you want the perturbed wavefunction
    psim = Psi(l1,l2)['psi']
    
    K = np.zeros((n,n))
    for i in range (n):
        for j in range(n): 
            if i != j: 
              K[i,j] = psi1[i,j]*psim[i,j]*(1/(np.abs(i-j)))
            else: 
              K[i,j] = 0 
            
    knum = np.trapz(np.trapz(K,x),y)
    
    kdenom = Psi(1,3)['energy'] - Psi(l1,l2)['energy']
    
    k = knum/kdenom
    
    return k 


#calculating perturbed wavefunction 
    
Psi_perturbed = K(1,2)*Psi(1,2)['psi'] + K(2,3)*Psi(2,3)['psi'] + K(1,4)*Psi(1,4)['psi'] + K(2,4)*Psi(2,4)['psi'] + K(3,4)*Psi(3,4)['psi'] 

#calculating perturbed prob density 
prob_density = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        prob_density[i,j] = np.abs(Psi_perturbed[i,j])**2

#well width 
W = 10*(10**-9) 

#step size  
n = 1000
#x-axis points and y-axis points 
x = np.linspace(0,W,n)
y = np.linspace(0,W,n)

#mesh points for plotting in 3-D 
xx, yy = np.meshgrid(x,y)
fig= plt.figure()
ax = fig.add_subplot(1, 2, 1, projection='3d')
p = ax.plot_surface(xx,yy,Psi_perturbed + Psi(1,3)['psi'], cmap='plasma')
ax.set_xlabel("$x_1$",fontsize = '20')
ax.set_ylabel("$x_2$",fontsize = '20')
ax.set_title('Wavefunction', fontsize = '20')

ax = fig.add_subplot(1, 2, 2, projection='3d')
p = ax.plot_surface(xx,yy,prob_density + Psi(1,3)['psi2'], cmap='plasma')
ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$")
ax.set_title('Probability Density')
