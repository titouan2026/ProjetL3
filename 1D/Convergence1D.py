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
L=700*10**(-8)

n0=10 #noeuds par longueur d'onde (100 pour avoir pas de soucis)

epsm=1 #coef mu et epsilon
epsd=3
mud=1
nm=1
nd=np.sqrt(epsd)
lambdad=lambda0/nd
rab=(nm-nd)/(nm+nd)
rbc=(nd-nm)/(nm+nd)
tab=2*nm/(nm+nd)
tbc=2*nd/(nm+nd)

k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))


#Initialisation de la géométrie de l'objet et du maillage

m1=0.1*L #debut de l'objet
m2=0.5*L #fin de l'objet
kd=k0*nd

rg=(rab+rbc*np.exp(2*1j*kd*(m2-m1)))/\
(1+rab*rbc*np.exp(2*1j*kd*(m2-m1))) #rapport réfléchie théorique
tg=(tab*tbc*np.exp(1j*kd*(m2-m1)))/\
(1+rab*rbc*np.exp(2*1j*kd*(m2-m1))) #rapport transmit théorique


R0=np.array([])
T0=np.array([])
NX=np.array([])

for i in range(2,21,1):
    n=np.round(n0*i**1.5)
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

    def S(x):
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

    ui=fi(maillage)

    for i in range(N+1):
        ud+=zeta[i]*Node(maillage,i)

    r = ud[ab]/ui[ab]
    t = (ud[bc]+ui[bc])/ui[bc]
    R0=np.append(R0,r)
    T0=np.append(T0,t)
    NX=np.append(NX,N+1)
    
R0-=rg
T0-=tg


plt.plot(np.log(NX),np.log(np.abs(T0)),color="red",linewidth=1,label="δ=t_calc-t_theo")
plt.plot(np.log(NX),np.log(np.abs(R0)),color="blue",linewidth=1,label="δ=r_calc-r_theo")
plt.scatter(np.log(NX),np.log(np.abs(R0)),s=5,marker="x",c="black")
plt.scatter(np.log(NX),np.log(np.abs(T0)),s=5,marker="x",c="black") 
plt.title("convergence 1D des coefficients de reflexion et transmition")
plt.xlabel("logN")
plt.ylabel("log|δ|")
plt.legend()
plt.show()


