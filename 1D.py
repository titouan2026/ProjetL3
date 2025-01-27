####libraries####

import numpy as np
import matplotlib.pyplot as plt

####initialisation variables####

lambda0=600*10**(-9)
theta=0
a=1
N=100
L=5*10**(-6)
h=L/N


epsm=1 #coef mu et epsilon
epsd=3+1j
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

def fiPrime(X):
    return a*k0*np.exp(X*k0*1j)*1j

def fiseconde(X):
    return -a*k0**(2)*np.exp(X*k0*1j)

ui=fi(maillage)


#Calcul du epsilon r et mu r

# def EpsilonContinue(x):
#     if x>= 50*h and x <= 65*h:
#         return epsd
#     else: return epsm
    

# def MuContinue(x):
#     if x>= 50*h and x <= 65*h:
#         return mud
#     else: return 1

def epsilon(X):
    a=np.min(np.where(grid==1))
    b=np.max(np.where(grid==1))
    return np.where((X/h>=a) & (X/h<=b),epsd,epsm)
    
    

def mu(X):
    a=np.min(np.where(grid==1))
    b=np.max(np.where(grid==1))
    return np.where((X/h>=a) & (X/h<=b),mud,1)


#Fonction de base nodale

# def node(i,X):  #fonction discrette
#     return np.where(X>=i,1,0)

def Node(X,i):  #fonction continue
    X=np.where((X>=h*(i-1)) & (X<=h*(i+1)),X,0)
    X=np.where((X>=h*(i-1)) & (X<=h*i),X/h-i+1,X)
    return np.where((X>=h*i) & (X<=h*(i+1)),-X/h+i+1,X)

def NodePrime(X,i):  #fonction continue
    X=np.where((X>=h*(i-1)) & (X<=h*(i+1)),X,0)
    X=np.where((X>=h*(i-1)) & (X<=h*i),1/h,X)
    return np.where((X>=h*i) & (X<=h*(i+1)),-1/h,X)

#Calcule de la fonction S

def Sapprox(x,i):
    return 1/mu(x)*fiPrime(x)*NodePrime(x,i)+k0**2*epsilon(x)*fi(x)*Node(x,i)

def S(x):
    return (1-1/mu(x))*fiseconde(x)+k0**2*(epsm-epsilon(x))*fi(x)


# def v(x,i):     #fonction 'continue'
#     if x >=h*(i-1) and x <= h*i:
#         return (x/h)-i+1
#     elif x >h*i and x <= h*(i+1):
#         return -(x/h)+i+1
#     else: return 0
    
# def vprime(x,i):
#     if x >=h*(i-1) and x <= h*i:
#         return 1/h
#     elif x >h*i and x <= h*(i+1):
#         return -1/h
#     else: return 0


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
    n=99+(N-2)*00
    dx=L/n
    x=np.arange(0,L,dx)
    f= NodePrime(x,i)*NodePrime(x,j)/mu(x) + k0**2*epsilon(x)*Node(x,i)*Node(x,j)
    return np.sum(f)*dx
    # for m in range(n):
    #     x=h*i+m*dx
    #     I+= (vprime(x,i)*vprime(x,j)/MuContinue(x) + k0**2*EpsilonContinue(x)*v(x,i)*v(x,j))*dx
    # return I

def IntegrateB(i):
    n=199
    dx=2*h/n
    if i == 0 :
        x=np.arange(0,h,dx)
        f=Sapprox(x,i)*Node(x,i)
        return np.sum(f)*dx +(fi(L)-fi(0))*k0*1j
    if i == N-1:
        x=np.arange(h*(N-1),h*N,dx)
        f=Sapprox(x,i)*Node(x,i)
        return np.sum(f)*dx +(fi(L)-fi(0))*k0*1j
        
        
    else:
        x=np.arange(h*(i-1),h*(i+1),dx)
        f=Sapprox(x,i)*Node(x,i)
        return np.sum(f)*dx 
        
    # for m in range(n):
    #     x=h*i+m*dx
    #     I+= S(x)*v(x,i)*dx

    # return I

B=np.array([])
A=np.zeros((N,N),dtype=complex)

for i in range(N):
    B=np.append(B,IntegrateB(i))
    
for i in range(N):
    for j in range(N):
        I=IntegrateA(i,j)
        A[j,i]=I
        A[i,j]=I

zeta = np.linalg.solve(A, B)

ud=0

n=99+(N-1)*100
dx=h*N/n
x=np.arange(0,L,dx)
ui=fi(x)
for i in range(N):
    ud+=zeta[i]*Node(x,i)
plt.plot(x,np.real(ud))
plt.plot(x,np.real(ui))
plt.plot([55*h,66*h],[0,0],linewidth=4,color="black")
plt.show()
plt.plot(x,np.real(ui+ud))
plt.show()
    
        
        
    