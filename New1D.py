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
L=0.4*10**(-5)

n=100 #noeuds par longueur d'onde (100 pour avoir pas de soucis)

epsm=1 #coef mu et epsilon
epsd=3+0.1j
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



#Champ initial

def fi(X: np.ndarray | float) -> np.ndarray | float:
    return a*np.exp(k0*X*1j)

def fiPrime(X: np.ndarray | float) -> np.ndarray | float:
    return a*k0*np.exp(X*k0*1j)*1j

def fiseconde(X: np.ndarray | float) -> np.ndarray | float:
    return -a*k0**2*np.exp(X*k0*1j)


#Epsilon et mu

def Epsilon(x: float) -> complex:
    
    if x >= m1 and x <= m2:
        return epsd
    else:
        return epsm

def Mu(x: float) -> complex:
    
    if x>= m1 and x <= m2:
        return mud
    else:
        return 1
    
def epsilon(X):
    
    return np.where((X>=m1) & (X<=m2),epsd,epsm)

def mu(X):
    
    return np.where((X>=m1) & (X<=m2),mud,1)


#Fonction de base nodale

def Node(X,i):
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
                
def NodePrime(X,i):  
    try: 
        if i == 0: raise IndexError
        X=np.where((X>=maillage[i-1]) & (X<=maillage[i+1]),X,0)
        X= np.where((X>maillage[i-1]) & (X< maillage[i]),1/(maillage[i]-maillage[i-1]),X)
        return np.where((X>maillage[i]) & (X<=maillage[i+1]),-1/(maillage[i+1]-maillage[i]),X)
    except:
        try:
            return np.where(X<=maillage[i+1],-1/(maillage[i+1]-maillage[i]),0)
        
        except:
            return np.where(X>maillage[i-1],1/(maillage[i]-maillage[i-1]),0)



#Calcule de la fonction S

def S(x: float) -> complex:
    return (1-1/Mu(x))*fiseconde(x)+k0**2*(epsm-Epsilon(x))*fi(x)



### création des matrices ###

K=np.zeros((N+1,N+1),dtype=complex)

for i in range(N):

    h=maillage[i+1]-maillage[i] #pas de maillage
    Khat=1/h*np.array([[1,-1],[-1,1]]) #matrice de raideur élémentaire
    K[i:i+2,i:i+2]+=Khat #assemblage de la matrice de raideur

M=np.zeros((N+1,N+1),dtype=complex)

for i in range(N):
    h=maillage[i+1]-maillage[i] #pas de maillage
    Mhat=h/6*np.array([[2,1],[1,2]]) #matrice de masse élémentaire
    M[i:i+2,i:i+2]+=k0**2*epsilon(maillage[i])*Mhat #assemblage de la matrice de masse

B=np.zeros(N+1,dtype=complex)

for i in range(N):
    h=maillage[i+1]-maillage[i] #pas de maillage
    Bhat=h/2*np.array([S(maillage[i]),S(maillage[i+1])]) #matrice de vecteur de charge élémentaire

    B[i:i+2]+=Bhat #assemblage du vecteur de charge
    
R=np.zeros((N+1,N+1),dtype=complex) #matrice de Robin

R[0,0]=k0*1j
R[N,N]=k0*1j

T = M-K+R #matrice finale


#Résolution du système et calcule du champ ud

zeta = np.linalg.solve(T, B)
ud=0

n=100*N
dx=L/n
x=np.arange(0,L+dx,dx)
ui=fi(x)

for i in range(N+1):
    ud+=zeta[i]*Node(x,i)


### Affichage du champ ###
    
plt.plot(x,np.real(ud),linewidth=1.2,label="ud",color="#5A9BD5")
plt.plot(x,np.real(ui),linewidth=1.2,label="ui",color="#E57373")

square_x = [m1, m1, m2, m2]
c=np.max(np.maximum(np.abs(np.real(ud)),np.abs(np.real(ui))))
square_y = [-c, c, c, -c] 
plt.fill(square_x, square_y, color="grey", alpha=0.2,label="Objet")
plt.text(m1+(m2-m1) / 20, 0.95*c , f'ε = {epsd}', horizontalalignment='left',
         verticalalignment='top', fontsize=7, color='black',weight='bold',alpha=0.8)    

plt.title("Champ d'une source infinie et d'une source objet")
plt.xlabel('x in metter')
plt.ylabel('Intensity')
plt.legend(loc='upper left')
plt.xlim(0, L)
plt.tight_layout()
plt.show()



plt.plot(x,np.real(ui+ud),linewidth=1.2,label="ud + ui",color="#5A9BD5")

c=np.max(ud+ui)
square_x = [m1, m1, m2, m2]
square_y = [-c, c, c, -c] 
plt.fill(square_x, square_y, color="grey", alpha=0.2,label="Objet")
plt.text(m1+(m2-m1) / 20, 0.95*c , f'ε = {epsd}', horizontalalignment='left',
         verticalalignment='top', fontsize=7, color='black',weight='bold',alpha=0.8)    

plt.title("Champ total")
plt.xlabel('x in metter')
plt.ylabel('Intensity')
plt.legend(loc='upper left')
plt.xlim(0, L)
plt.tight_layout()
plt.show()
