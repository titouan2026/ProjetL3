####libraries####

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from matplotlib import cm
import math
import mesh_io as msh
import meshio

### DebugGing tool ###

# import sys

# if not sys.warnoptions:  #Raise warnings as error
#     import warnings
#     warnings.simplefilter("error")

####initialisation variables####
    
lambda0=500*10**(-9)
theta=np.pi/4
a=1
R=1*10**(-6)

n=1 #noeuds par longueur d'onde (>50 pour avoir pas de soucis)

epsm=1 #coef mu et epsilon
epsd=3
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
    mesh=msh.GetData(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4)
    print("mesh found")
    
except:
    print("creating mesh")
    msh.CreateMesh(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4)
    mesh=msh.GetData(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4)


# Charger le maillage Gmsh


# Extraire les coordonnées des points
points = mesh.points  # Liste des coordonnées (x, y) des nœuds
cells = mesh.cells_dict["triangle"]

def shape_functions(P1, P2, P3):
    """ Calcule les coefficients des fonctions nodales P1 sur un triangle. """
    x1, y1 = P1
    x2, y2 = P2
    x3, y3 = P3

    # Calcul de l'aire du triangle
    A = 0.5 * abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))

    # Coefficients des fonctions de forme
    phi1 = [(x2*y3 - x3*y2) / (2*A), (y2 - y3) / (2*A), (x3 - x2) / (2*A)]
    phi2 = [(x3*y1 - x1*y3) / (2*A), (y3 - y1) / (2*A), (x1 - x3) / (2*A)]
    phi3 = [(x1*y2 - x2*y1) / (2*A), (y1 - y2) / (2*A), (x2 - x1) / (2*A)]
    
    return phi1, phi2, phi3

# Test sur le premier élément du maillage
print(points[cells[0][0]][:2], points[cells[0][1]][:2], points[cells[0][2]][:2])
P1, P2, P3 = points[cells[0][0]][:2], points[cells[0][1]][:2], points[cells[0][2]][:2]
phi1, phi2, phi3 = shape_functions(P1, P2, P3)

print("φ1(x, y) = ", phi1)
print("φ2(x, y) = ", phi2)
print("φ3(x, y) = ", phi3)


# x=np.arange(-R,R,2*R/1000)
# y=np.arange(-R,R,2*R/1000)
# xgrid, ygrid = np.meshgrid(x, y)
# Z=Node2D(xgrid,ygrid,31)

# fig, ax = plt.subplots()
# mask = (xgrid**2 + ygrid**2) <= (R**2)
# Z = np.where(mask, Z, np.nan)

# contour = ax.contourf(xgrid, ygrid, Z, cmap=cm.hot)

# plt.scatter(X,Y, color='black', s=1)
# plt.colorbar(contour, label='Z')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Contour plot with color mapping')

# # Afficher le graphique
# plt.show()


## donées aléatoires ##
# nx = 10*int(N**0.4) #peut créé des problèmes si N est trop grand
# xg = np.linspace(X.min(), X.max(), nx)
# yg = np.linspace(Y.min(), Y.max(), nx)
# xgrid, ygrid = np.meshgrid(xg, yg)
# ctr_f = griddata((X, Y), Z, (xgrid, ygrid), method='linear')

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1) 
# ax.contourf(xgrid, ygrid, ctr_f, cmap=cm.autumn)

# circle = plt.Circle((0, 0), R*0.985, color='black', fill=False)
# plt.scatter(X,Y, color='black', s=1)
# ax.set_aspect('equal', adjustable='box')
# ax.set_xlim(-R, R)
# ax.set_ylim(-R, R)

# # Ajout du cercle aux axes
# ax.add_artist(circle)

# plt.show()

