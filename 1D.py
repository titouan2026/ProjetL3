####libraries####

import numpy as np
import matplotlib.pyplot as plt

####initialisation variables####

lambda0=525*10**(-9)
theta=0
a=1
N=100
L=5*10**(-6)
h=L/N


epsm=1 #coef mu et epsilon
epsd=80
mud=1

k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))


#Initialisation de la géométrie de l'objet et du maillage

maillage = np.arange(N)*h
grid = np.zeros(N)
grid[50:66]=1



#Calcul du champ initial

def fi(X):
    return a*np.exp(k0*X*1j)

def fiseconde(X):
    return -a*k0**(2)*np.exp(X*k0*1j)

ui=fi(maillage)


#Calcule du epsilon r et mu r

def EpsilonContinue(x):
    if x>= 50*h and x <= 65*h:
        return epsd
    else: return epsm
    

def MuContinue(x):
    if x>= 50*h and x <= 65*h:
        return mud
    else: return 1

def epsilon():
    
    return epsd*grid+(1-grid)*epsm
    

def mu():
    return mud*grid+(1-grid)


#Calcule de la fonction S

def S(x):
    return (1-1/MuContinue(x))*fiseconde(x)+k0**2*(epsm-EpsilonContinue(x))*fi(x)


#Fonction de base nodale

def node(i,X):  #fonction discrette
    return np.where(X==i,1,0)

def v(x,i):     #fonction 'continue'
    if x >=h*(i-1) and x <= h*i:
        return (x/h)-i+1
    elif x >h*i and x <= h*(i+1):
        return -(x/h)+i+1
    else: return 0
    
def vprime(x,i):
    if x >=h*(i-1) and x <= h*i:
        return 1/h
    elif x >h*i and x <= h*(i+1):
        return -1/h
    else: return 0


### Intégration ###

def IntegrateA(i,j):
    d=np.abs(i-j)
    
    if i>j: #la matrice est symétrique on traitera donc ce cas apres
        return 0 
    else:
        pass
    
    if d>=2:
        return 0
    else:
        pass
    
    n=100*(d+1)
    dx=(h*(d+2))/n
    I=0
    for m in range(n):
        x=h*i+m*dx
        I+= (vprime(x,i)*vprime(x,j)/MuContinue(x) + k0**2*EpsilonContinue(x)*v(x,i)*v(x,j))*dx
    return I

def IntegrateB(i):
    n=100
    dx=2*h/n
    I=0
    for m in range(n):
        x=h*i+m*dx
        I+= S(x)*v(x,i)*dx

    return I

B=np.array([])
A=np.zeros((N-1,N-1))

for i in range(1,N):
    B=np.append(B,IntegrateB(i))
    
for i in range(1,N):
    for j in range(i,N):
        I=IntegrateA(i,j)
        A[j-1,i-1]=I
        A[i-1,j-1]=I

zeta = np.linalg.solve(A, B)

ud=0
for i in range(1,N):
    ud+=zeta[i-1]*node(i,np.arange(N))

plt.plot(np.arange(N),ud)
plt.plot(np.arange(N),ui)
plt.show()
plt.plot(np.arange(N),ui+ud)
plt.show()
    
        
        
    