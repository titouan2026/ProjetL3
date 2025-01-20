####libraries####

import numpy as np
import matplotlib.pyplot as plt

####initialisation variables####

lambda0=550*10**(-9)
theta=0
A=1
N=21
L=5*10**(-7)
h=L/N


epsm=1 #coef mu et epsilon
epsd=80
mud=1

k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))




#initialisation de la géométrie de l'objet

grid = np.zeros(21)
grid[5:16]=1

#calcul du champ initial

def fi(x):
    return A*np.exp(k0*x*1j)

def fiseconde(x):
    return -A*k0**(2)*np.exp(x*k0*1j)

ui=fi(np.arange(21)*h)

#calcule du epsilon r et mu r

def epsilon(x):
    
    return epsd*grid+(1-grid)*epsm
    

def mu(x):
    return mud*grid+(1-grid)
    
def S(x):
    return (1-1/mu(x))*fiseconde(x)+k0**2*(epsm-epsilon(x))*fi(x)

s=S(np.arange(N)*h)

def base(i,x):
    return np.where(x==i,1,0)
