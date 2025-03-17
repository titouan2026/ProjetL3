import gmsh_api.gmsh as gmsh
import numpy as np

def CreateMesh(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4):
    
    lc=1
    a=0.2*R
    b=-0.2*R
    c=np.real(nd)/nm
    gmsh.finalize()
    gmsh.initialize()
    
    gmsh.option.setNumber("General.Verbosity", 2)

    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lambda0/n)

    circle=gmsh.model.occ.addCircle(0.0, 0.0, 0.0,R,5,angle1=0.,angle2=2*np.pi)
    
    point1 = gmsh.model.occ.add_point(x1, y1, 0, lc/n)
    point2 = gmsh.model.occ.add_point(x2, y2, 0, lc/n)
    point3 = gmsh.model.occ.add_point(x3, y3, 0, lc/n)
    point4 = gmsh.model.occ.add_point(x4, y4, 0, lc/n)

    line1 = gmsh.model.occ.add_line(point1, point2)
    line2 = gmsh.model.occ.add_line(point2, point3)
    line3 = gmsh.model.occ.add_line(point3, point4)
    line4 = gmsh.model.occ.add_line(point4, point1)


    gmsh.model.occ.addCurveLoop([line1,line2,line3,line4] ,7)

    gmsh.model.occ.addCurveLoop([circle] ,6)

    gmsh.model.occ.addPlaneSurface([6, 7],8)

    gmsh.model.occ.addPlaneSurface([7],9) #Suface du millieu, peut etre modifié mais dois resté numéo 9!


    gmsh.model.occ.synchronize()
    field_size = gmsh.model.mesh.field.add("MathEval",13)
    gmsh.model.mesh.field.setString(field_size, "F", "{}".format(lambda0/n))

    field_size2 = gmsh.model.mesh.field.add("MathEval",14)
    gmsh.model.mesh.field.setString(field_size2, "F", "{}".format(lambda0/(n*c)))

    field_restrict = gmsh.model.mesh.field.add("Restrict",15)
    gmsh.model.mesh.field.setNumbers(field_restrict, "SurfacesList", [9])
    gmsh.model.mesh.field.setNumber(field_restrict, "IField", field_size2)

    gmsh.model.mesh.field.add("Min", 12)
    gmsh.model.mesh.field.setNumbers(12, "FieldsList",[15,13])
    gmsh.model.mesh.field.setAsBackgroundMesh(12)

    gmsh.model.mesh.generate(2)

    gmsh.write("C:/Users/PC/Desktop/cours/ProjetL3/Meshes/{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.msh".format(R, nm, nd, lambda0, n, x1, y1, x2, y2, x3, y3, x4, y4))

    gmsh.finalize()

def GetData(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4):
    
    dim=-1
    tag=-1
    gmsh.initialize()
    
    gmsh.option.setNumber("General.Verbosity", 2)
    
    gmsh.open("C:/Users/PC/Desktop/cours/ProjetL3/Meshes/{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}.msh".format(R,nm,nd,lambda0,n,x1,y1,x2,y2,x3,y3,x4,y4))
    nodeTags, coords, parametricCoord= gmsh.model.mesh.getNodes(dim,tag)
    coords= coords.reshape((-1,3))[:,0:2]
    gmsh.finalize()
    return coords
