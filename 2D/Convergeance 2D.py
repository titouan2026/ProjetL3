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
import mesh_ioconv as msh
import meshio
import permitivity as perm 

### DebugGing tool ###

# import sys

# if not sys.warnoptions:  #Raise warnings as error
#     import warnings
#     warnings.simplefilter("error")

####initialisation variables####
    
lambda0=1000e-9
theta=0
a=1
L=500e-8
l=100e-9



#noeuds par longueur d'onde (100 pour avoir pas de soucis)

epsm=1 #coef mu et epsilon
epsd=3
mud=1
nm=1
nd=np.sqrt(epsd)
lambdad=lambda0/nd
lambdadreel=lambda0/np.real(nd)
rab=(nm-nd)/(nm+nd)
rbc=(nd-nm)/(nm+nd)
tab=2*nm/(nm+nd)
tbc=2*nd/(nm+nd)

k0=2*np.pi/lambda0
k=(k0*np.cos(theta),k0*np.sin(theta))
kd=k0*nd


m1=0.1*L #debut de l'objet
m2=0.4*L #fin de l'objet
kd=k0*nd

rg=(rab+rbc*np.exp(2*1j*kd*(m2-m1)))/\
(1+rab*rbc*np.exp(2*1j*kd*(m2-m1))) #rapport réfléchie théorique
tg=(tab*tbc*np.exp(1j*kd*(m2-m1)))/\
(1+rab*rbc*np.exp(2*1j*kd*(m2-m1))) #rapport transmit théorique
rg=-1.0861267821657767+0.0019820742306421594j
tg=-0.08610269083951483+0.001980641291361989j

n0=2

Ref=np.array([])
Trans=np.array([])
NX=np.array([])

for g in range(3,11,2):
    L=400e-8
    l=100e-9
    n=np.round(n0*g**1.5) 
    try:
        maillage=msh.GetData(L,l,nm,nd,lambda0,n)
        print("mesh found")
        
    except:
        print("creating mesh")
        msh.CreateMesh(L,l,nm,nd,lambda0,n)
        maillage=msh.GetData(L,l,nm,nd,lambda0,n)
        print("mesh created")


    # Extraire les coordonnées des points
    points = maillage.points[:,:2]  # Liste des coordonnées (x, y) des nœuds
    cells = np.array(maillage.cells_dict["triangle"])
    N=np.size(points)

    X=points[:,0]
    Y=points[:,1]


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
                
    
    Rob=np.zeros((length,length),dtype=complex)


    for line in boundaries_group:
        x1,y1=points[line[0]]
        x2,y2=points[line[1]]
        
        loctoglb=[line[0],line[1]]
        for i in range(2):
            for j in range(2):
                if Rob[loctoglb[i],loctoglb[j]]==0: 
                    Rob[loctoglb[i],loctoglb[j]]+=k0*1j
                else: pass
                
        

    T=Rob-A1+M
    zeta = spsolve(csc_matrix(T), B)
    # zeta = spsolve(csc_matrix(T), B)
    Z=np.zeros(length,dtype=complex)

    for i in range(length):

        Z[i]=zeta[i]
    ui= fi(points.T)
    
    # triang = mtri.Triangulation(X, Y)
    # # # Tracé avec plt.tricontourf
    # fig, ax = plt.subplots()
    # contour = ax.tricontourf(triang, np.real(ui), cmap=cm.twilight)
    # plt.colorbar(contour, label='Z')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.title('Triangulated contour plot')
    # plt.show()
    
    
    nx = 2000 #peut créé des problèmes si nx est trop grand
    ny=nx
    dx= (X.max()-X.min())/nx
    

    xg = np.linspace(X.min(), X.max(), nx)
    yg = np.linspace(Y.min(), Y.max(), ny) 
    xgrid, ygrid = np.meshgrid(xg, yg)
    # Interpolation des données sur la grille
    ud = griddata((X,Y), Z, (xgrid, ygrid), method='linear')
    ui= griddata((X,Y), ui, (xgrid, ygrid), method='linear')
    
    
    ab=round(0.2*L/dx)
    bc=round(0.6*L/dx)
    
    # r = ud[1000][ab]/ui[1000][ab]
    # t = (ud[1000][bc]+ui[1000][bc])/ui[1000][bc]
    
    r = ud[ab][1000]/ui[ab][1000]
    t = (ud[bc][1000]+ui[bc][1000])/ui[bc][1000]
    
    Ref=np.append(Ref,r)
    Trans=np.append(Trans,t)
    NX=np.append(NX,N+1)
    
    print(r,t,rg,tg) 
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1) 
    # ax.contourf(xgrid, ygrid, np.real(ui+ud), cmap=cm.autumn)

    # ax.set_aspect('equal', adjustable='box')


    # plt.show()
Ref-=rg
Trans-=tg

plt.plot(np.log(NX),np.log(np.abs(Ref)),color="red",linewidth=1,label="δ=t_calc-t_theo")
plt.plot(np.log(NX),np.log(np.abs(Trans)),color="blue",linewidth=1,label="δ=r_calc-r_theo")
plt.scatter(np.log(NX),np.log(np.abs(Trans)),s=5,marker="x",c="black")
plt.scatter(np.log(NX),np.log(np.abs(Ref)),s=5,marker="x",c="black") 
plt.title("convergence 2D des coefficients de reflexion et transmition")
plt.xlabel("logN")
plt.ylabel("log|δ|")
plt.legend()
plt.show()

    
    

