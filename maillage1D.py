####libraries####

import numpy as np
import matplotlib.pyplot as plt
import math

### DebugGing tool ###

# import sys

# if not sys.warnoptions:  #Raise warnings as error
#     import warnings
#     warnings.simplefilter("error")

####initialisation variables####
    
lambda0=700*10**(-9)
theta=0
a=1
L=1.8*10**(-6)

n=5 #noeuds par longueur d'onde (100 pour avoir pas de soucis)

epsm=1 #coef mu et epsilon
epsd=3+1j
mud=1
nm=1
nd=np.sqrt(epsd)
lambdad=lambda0/nd


k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))
kd=k0*nd


#Initialisation de la géométrie de l'objet et du maillage

m1=0.3*L
m2=0.7*L

maillage = np.concatenate((np.linspace(0, m1, math.ceil(m1/(lambda0/n)),endpoint=False), np.linspace(m1, m2, math.ceil((m2-m1)/(lambdad/n)),endpoint=False), np.linspace(m2, L, math.ceil((L-m2)/(lambda0/n)),endpoint=True)))


ab=np.where(maillage==m1)[0][0] #indice (python) des interfaces
bc=np.where(maillage==m2)[0][0]
N=np.size(maillage)-1
grid = np.zeros(N+1)
grid[ab:bc+1]=1



#Calcul du champ initial

def fi(X):
    return a*np.exp(k0*X*1j)

def fiPrime(X):
    return a*k0*np.exp(X*k0*1j)*1j

def fiseconde(X):
    return -a*k0**2*np.exp(X*k0*1j)



def epsilon(X):
    
    return np.where((X>=m1) & (X<=m2),epsd,epsm)

def mu(X):
    
    return np.where((X>=m1) & (X<=m2),mud,1)


def Epsilon(X):
    
    if X >= m1 and X <= m2:
        return epsd
    else:
        return epsm

def Mu(X):
    
    if X>= m1 and X <= m2:
        return mud
    else:
        return 1



#Fonction de base nodale

def Node(X,i):  #fonction continue
    try:
        if i == 0: raise IndexError
        X=np.where((X>=maillage[i-1]) & (X<=maillage[i+1]),X,0)
        X= np.where((X>=maillage[i-1]) & (X<= maillage[i]),(X-maillage[i-1])/(maillage[i]-maillage[i-1]),X)
        return np.where((X>maillage[i]) & (X<=maillage[i+1]),(maillage[i+1]-X)/(maillage[i+1]-maillage[i]),X)
    except:
        try:
            return np.where(X<=maillage[i+1],(maillage[i+1]-X)/(maillage[i+1]-maillage[i]),0)
        except:
            return np.where(X>=maillage[i-1] ,(X-maillage[i-1])/(maillage[i]-maillage[i-1]),0)
            
        
def NodePrime(X,i):  #fonction continue
     
    try: 
        if i == 0: raise IndexError
        X=np.where((X>=maillage[i-1]) & (X<=maillage[i+1]),X,0)
        X= np.where((X>maillage[i-1]) & (X< maillage[i]),1/(maillage[i]-maillage[i-1]),X)
        return np.where((X>maillage[i]) & (X<=maillage[i+1]),-1/(maillage[i+1]-maillage[i]),X)
    
    except:
        print('cc')
        try:
            return np.where(X<=maillage[i+1],-1/(maillage[i+1]-maillage[i]),0)
        
        except:
            print('cc')
            return np.where(X>maillage[i-1],1/(maillage[i]-maillage[i-1]),0)



#Calcule de la fonction S

def S(x):
    return (1-1/Mu(x))*fiseconde(x)+k0**2*(epsm-Epsilon(x))*fi(x)





A=np.where(grid==1)[0][0]
B=np.where(grid==1)[0][-1]

plt.plot(maillage,Node(maillage,10),linewidth=1,color='black',label='10ième fonction nodale')
plt.scatter(maillage,[0]*np.size(maillage),color='black',s=10,label='Noeuds du millieux')
plt.scatter(maillage[A:B+1],[0]*np.size(np.where(grid==1)),color='red',s=10,label="Noeuds de l'objet")
plt.ylim(-0.2,1.6)
c=1.5
square_x = [m1, m1, m2, m2]
square_y = [-c, c, c, -c] 
plt.fill(square_x, square_y, color="grey", alpha=0.2,label="Objet")
plt.text(m2-(m2-m1) / 20, 0.95*c , f'ε = {epsd}', horizontalalignment='right', verticalalignment='top', fontsize=7, color='black',weight='bold',alpha=0.8)    
plt.xlabel('x in metter')
plt.ylabel('f(x)')
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()