####libraries####

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from matplotlib import cm
import math
import mesh_io as msh
import meshio
import permitivity as perm 

### DebugGing tool ###

# import sys

# if not sys.warnoptions:  #Raise warnings as error
#     import warnings
#     warnings.simplefilter("error")

#######################

####initialisation variables####
    
lambda0=500e-9
theta=np.pi/6
a=1
R=1.1e-6

n=10 #noeuds par longueur d'onde (100 pour avoir pas de soucis)

epsm=1 #coef mu et epsilon
epsd=9+1j
mud=1
nm=1
nd=np.sqrt(epsd)
lambdad=lambda0/nd
lambdadreel=lambda0/np.real(nd)

k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))
kd=k0*nd


#Initialisation de la géométrie de l'objet et du maillage

#coordonnées des points de l'objet 
x1=0.3*R
y1=0.2*R
x2=0.3*R
y2=-0.2*R
x3=-0.3*R
y3=-0.2*R
x4=-0.3*R
y4=0.2*R

try:
    maillage=msh.GetData(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4)
    print("mesh found")
    
except:
    print("creating mesh")
    msh.CreateMesh(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4)
    maillage=msh.GetData(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4)
    print("mesh created")


# Extraire les coordonnées des points
points = maillage.points[:,:2]  # Liste des coordonnées (x, y) des nœuds
cells = maillage.cells_dict["triangle"]

X=points[:,0]
Y=points[:,1]
# XY=np.meshgrid(X,Y)

# Sortedx=np.sort(np.unique(X))
# Sortedy=np.sort(np.unique(Y))
# N=len(X)

# Liste pour stocker les nœuds associés à ce groupe
object_group = []    
boundaries_group = []

# Parcourir les cellules (éléments maillages)
physical_tags = maillage.get_cell_data("gmsh:physical", "triangle")

for i,tag in enumerate(physical_tags):
        
    if tag==17:
        object_group.append(cells[i])
    else:
        pass

line_cells = maillage.cells_dict["line"]
physical_tags = maillage.get_cell_data("gmsh:physical", "line")


for i, tag in enumerate(physical_tags):
    
    if tag==16:
        boundaries_group.append(line_cells[i])
    else:
        pass
        
            
object_group=np.array(object_group)
boundaries_group=np.array(boundaries_group)



#Champ initial

def fi(XY: np.ndarray) -> np.ndarray | float:
    return a*np.exp(np.dot(k,XY)*1j)

def fiPrime(XY: np.ndarray) -> np.ndarray | float:
    return a*k0*np.exp(np.dot(k,XY)*1j)*1j

def fiseconde(XY: np.ndarray) -> np.ndarray | float:
    return -a*k0**2*np.exp(np.dot(k,XY) *k0*1j)
            
def area(x1, y1, x2, y2, x3, y3):
    return abs((x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)) / 2.0)

def scale_coordinates(factor, *coords):
    return [coord * factor for coord in coords]

def isIn(cell):
    if cell in object_group:
        return True
    else:   
        return False

def epsilon(cell) -> complex:
    if isIn(cell):
        return epsd
    else:
        return epsm
    
def mu(cell) -> complex:
    if isIn(cell):
        return mud
    else:
        return 1

def S(x,y: float,cell) -> complex:
    return (1-1/mu(cell))*fiseconde((x,y))+k0**2*(epsm-epsilon(cell))*fi((x,y))


### Création des matrices et vecteurs ###

length=len(points)

M0=(1/12)*np.array([[2,1,1],[1,2,1],[1,1,2]]) 
M=np.zeros((length,length),dtype=complex)

for i in range(len(cells)):
    x1,y1=points[cells[i][0]]
    x2,y2=points[cells[i][1]]
    x3,y3=points[cells[i][2]]
    A=area(x1,y1,x2,y2,x3,y3)
    loctoglb=[cells[i][0],cells[i][1],cells[i][2]]
    tempEps=epsilon(cells[i])
    for i in range(3):
        for j in range(3):
            M[loctoglb[i],loctoglb[j]]+=A*k0**2*tempEps*M0[i,j]

B=np.zeros(length,dtype=complex)

for i in range(len(cells)):
    x1,y1=points[cells[i][0]]
    x2,y2=points[cells[i][1]]
    x3,y3=points[cells[i][2]]
    A=area(x1,y1,x2,y2,x3,y3)
    loctoglb=[cells[i][0],cells[i][1],cells[i][2]]
    for j in range(3):
        B[loctoglb[j]]+=A/3*S(points[loctoglb[j]][0],points[loctoglb[j]][1],cells[i])

A1=np.zeros((length,length),dtype=complex)

for i in range(len(cells)):
    x1,y1=points[cells[i][0]]
    x2,y2=points[cells[i][1]]
    x3,y3=points[cells[i][2]]
    A=area(x1,y1,x2,y2,x3,y3)
    # a=1/(2*A)*np.array([(x2*y3-x3*y2),(x3*y1-x1*y3),(x1*y2-x2*y1)])
    b=1/(2*A)*np.array([(y2-y3),(y3-y1),(y1-y2)])
    c=1/(2*A)*np.array([(x3-x2),(x1-x3),(x2-x1)])
    A1k=np.array([[b[0]**2+c[0]**2,b[0]*b[1]+c[0]*c[1],b[0]*b[2]+c[0]*c[2]],
                  [b[1]*b[0]+c[1]*c[0],b[1]**2+c[1]**2,b[1]*b[2]+c[1]*c[2]],
                  [b[2]*b[0]+c[2]*c[0],b[2]*b[1]+c[2]*c[1],b[2]**2+c[2]**2]])
    loctoglb=[cells[i][0],cells[i][1],cells[i][2]] 
    for i in range(3):
        for j in range(3):
            A1[loctoglb[i],loctoglb[j]]+=A*A1k[i,j]
            
      

R0=(1/6)*np.array([[2,1],[1,2]]) ### a revoir !!!!!!!!!!!!!!!!
Rob=np.zeros((length,length),dtype=complex)


for line in boundaries_group:
    x1,y1=points[line[0]]
    x2,y2=points[line[1]]
    L=math.sqrt((x2-x1)**2+(y2-y1)**2)
    loctoglb=[line[0],line[1]]
    for i in range(2):
        for j in range(2):
            Rob[loctoglb[i],loctoglb[j]]+=np.abs(L)*k0*1j*R0[i,j]
    

T=Rob-A1+M
zeta = spsolve(csc_matrix(T), B)
Z=np.zeros(length,dtype=complex)

for i in range(length):

    Z[i]=zeta[i]
Z+= fi(points.T)


print("Start plotting")
 
x1=0.3*R
y1=0.2*R
x2=0.3*R
y2=-0.2*R
x3=-0.3*R
y3=-0.2*R
x4=-0.3*R
y4=0.2*R

nx = 2000 #peut créé des problèmes si nx est trop grand
ny=nx

xg = np.linspace(X.min(), X.max(), nx)
yg = np.linspace(Y.min(), Y.max(), ny)
xgrid, ygrid = np.meshgrid(xg, yg)
# Interpolation des données sur la grille
zgrid = griddata((Y,X), np.real(Z), (xgrid, ygrid), method='linear')
# Tracé avec plt.imshow
fig =plt.figure()
ax = fig.add_subplot(1, 1, 1) 
plt.imshow(zgrid.T, extent=(X.min(), X.max(), Y.min(), Y.max()), origin='lower',vmin=np.real(Z).min()*1.2*a, vmax=np.real(Z).max()*1.1*a, cmap='hot')

# Ajouts
square_x = [x1, x2, x3, x4, x1]
square_y = [y1, y2, y3, y4, y1]
plt.plot(square_x, square_y, color='black')
ax.set_xlim(-R*1.1, R*1.1)
ax.set_ylim(-R*1.1, R*1.1)
circle = plt.Circle((0, 0), R*0.995, color='black', fill=False)
ax.add_artist(circle)
ax.set_aspect('equal', adjustable='box')

plt.colorbar(label='Z')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Image plot with interpolation')
plt.show()



### Autres méthodes de tracages (surement non fonctionelles) ###


#création d'une grille régulière pour l'affichage
# nx = 10*int(N**0.4) #peut créé des problèmes si N est trop grand
# xg = np.linspace(X.min(), X.max(), nx)
# yg = np.linspace(Y.min(), Y.max(), nx)
# xgrid, ygrid = np.meshgrid(xg, yg)
# ctr_f = griddata((X, Y), Z, (xgrid, ygrid), method='linear')

# mask = (xgrid**2 + ygrid**2) <= ((R*0.985)**2)
# ctr_f = np.where(mask, ctr_f, np.nan)

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1) 
# ax.contourf(xgrid, ygrid, ctr_f, cmap=cm.autumn)

# circle = plt.Circle((0, 0), R*0.985, color='black', fill=False)
# ax.set_aspect('equal', adjustable='box')
# ax.set_xlim(-R, R)
# ax.set_ylim(-R, R)

# # Ajout du cercle aux axes
# ax.add_artist(circle)

# plt.show()





# triang = mtri.Triangulation(X, Y)

# # Tracé avec plt.tricontourf
# fig, ax = plt.subplots()
# contour = ax.tricontourf(triang, np.real(Z), cmap=cm.twilight)
# plt.colorbar(contour, label='Z')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Triangulated contour plot')

# # Ajout du cercle aux axes
# circle = plt.Circle((0, 0), R*0.985, color='black', fill=False)
# ax.add_artist(circle)
# ax.set_aspect('equal', adjustable='box')
# ax.set_xlim(-R, R)
# ax.set_ylim(-R, R)

# plt.show()

