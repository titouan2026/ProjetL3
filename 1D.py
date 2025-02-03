####libraries####

import numpy as np
import matplotlib.pyplot as plt

####initialisation variables####

lambda0=700*10**(-9)
theta=0
a=1
N=300
L=5*10**(-6)
h=L/N


epsm=1 #coef mu et epsilon
epsd=3
mud=1

k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))


#Initialisation de la géométrie de l'objet et du maillage

maillage = np.arange(N+1)*h
grid = np.zeros(N+1)
grid[100:250]=1
m1=np.min(np.where(grid==1))
m2=np.max(np.where(grid==1))


#Calcul du champ initial

def fi(X):
    return a*np.exp(k0*X*1j)

def fiPrime(X):
    return a*k0*np.exp(X*k0*1j)*1j

def fiseconde(X):
    return -a*k0**2*np.exp(X*k0*1j)




#Calcul du epsilon r et mu r

def epsilon(X):
    
    return np.where((X/h>=m1) & (X/h<=m2),epsd,epsm)

    
    
def mu(X):
    
    return np.where((X/h>=m1) & (X/h<=m2),mud,1)


#Fonction de base nodale

def Node(X,i):  #fonction continue
    X=np.where((X>=h*(i-1)) & (X<=h*(i+1)),X,0)
    X= np.where((X>=h*(i-1)) & (X< h*i),(X-h*(i-1))/h,X)
    return np.where((X>h*i) & (X<=h*(i+1)),(h*(i+1)-X)/h,X)

def NodePrime(X,i):  #fonction continue
    X = np.where((X>=h*(i-1)) & (X<=h*(i+1)),X,0)
    X = np.where((X>h*(i-1)) & (X < (h*i)),1/h,X)
    return np.where((X>h*i) & (X<=h*(i+1)),-1/h,X)


#Calcule de la fonction S

def Sapprox(x,i):
    return 1/mu(x)*fiPrime(x)*NodePrime(x,i)+k0**2*epsilon(x)*fi(x)*Node(x,i)

def S(x):
    return (1-1/mu(x))*fiseconde(x)+k0**2*(epsm-epsilon(x))*fi(x)

### Intégration ###

def IntegrateA(i,j):
    d=np.abs(i-j)
    
    if d>=2:
        return 0
    else:
        pass
    n=99+(N-1)*100
    dx=L/n
    x=np.arange(0,L,dx,dtype=complex)
    f= NodePrime(x,i)*NodePrime(x,j)/mu(x) + k0**2*epsilon(x)*Node(x,i)*Node(x,j)
    return np.sum(f)*dx


def IntegrateB(i):
    n=199
    dx=2*h/n
    if i == 0 :
        dx=h/99
        x=np.arange(0,h,dx,dtype=complex)
        f=S(x)*Node(x,i)
        return np.sum(f,dtype=complex)*dx
    if i == N:
        dx=h/99
        x=np.arange(L-h,L,dx,dtype=complex)
        f=S(x)*Node(x,i)
        return np.sum(f,dtype=complex)*dx
        
        
    else:
        x=np.arange(h*(i-1),h*(i+1),dx)
        f=S(x)*Node(x,i)
        return np.sum(f,dtype=complex)*dx 
        

B=np.array([])
A=np.zeros((N+1,N+1),dtype=complex)

for i in range(N+1):
    B=np.append(B,IntegrateB(i))
    
for i in range(N+1):
    for j in range(N+1):
        I=IntegrateA(i,j)
        A[i,j]=I
A[0,0] += k0*1j
A[N,N] += k0*1j
zeta = np.linalg.solve(A, B)


ud=0

n=99+(N-1)*100
dx=L/n
x=np.arange(0,L+dx,dx)
ui=fi(x)

for i in range(N+1):
    ud+=zeta[i]*Node(x,i)
    
   
    
plt.plot(x,np.real(ud),linewidth=1.2,label="ud",color="#5A9BD5")
plt.plot(x,np.real(ui),linewidth=1.2,label="ui",color="#E57373")

square_x = [m1*h, m1*h, m2*h, m2*h]
c=np.max(np.maximum(ud,ui))
square_y = [-c, c, c, -c] 
plt.fill(square_x, square_y, color="grey", alpha=0.2,label="Objet")

plt.title("Champ d'une source infinie et d'une source objet")
plt.xlabel('x in metter')
plt.ylabel('Intensity')
plt.legend(loc='upper left')
plt.xlim(0, L)
plt.tight_layout()
plt.show()



plt.plot(x,np.real(ui+ud),linewidth=1.2,label="ud + ui",color="#5A9BD5")

square_x = [m1*h, m1*h, m2*h, m2*h]
square_y = [-np.max(ud+ui), np.max(ud+ui), np.max(ud+ui), -np.max(ud+ui)] 
plt.fill(square_x, square_y, color="grey", alpha=0.2,label="Objet")

plt.title("Champ total")
plt.xlabel('x in metter')
plt.ylabel('Intensity')
plt.legend(loc='upper left')
plt.xlim(0, L)
plt.tight_layout()
plt.show()


    
        
        
    