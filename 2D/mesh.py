import gmsh_api.gmsh as gmsh
import numpy as np
import shelve

def CreatingMesh(R,nm,nd,lambda0,n):

    lc=1
    a=0.2*R
    b=-0.2*R
    c=np.real(nd)/nm

    gmsh.initialize()

    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lambda0/n)

    circle=gmsh.model.occ.addCircle(0.0, 0.0, 0.0,R,5,angle1=0.,angle2=2*np.pi)
    point1 = gmsh.model.occ.add_point(1.5*a, a, 0, lc/n)
    point2 = gmsh.model.occ.add_point(1.5*a, b, 0, lc/n)
    point3 = gmsh.model.occ.add_point(1.5*b, b, 0, lc/n)
    point4 = gmsh.model.occ.add_point(1.5*b, a, 0, lc/n)

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
    gmsh.fltk.run()
    # gmsh.write("Meshes\Circle.msh")
    gmsh.write("Meshes/{}_{}_{}_{}_{}.msh".format(R,nm,nd,lambda0,n))
    gmsh.finalize()

def GetData(R,nm,nd,lambda0,n):
    dim=-1
    tag=-1
    gmsh.initialize()
    
    gmsh.open("Meshes/{}_{}_{}_{}_{}.msh".format(R,nm,nd,lambda0,n))
    nodeTags, coords, parametricCoord= gmsh.model.mesh.getNodes(dim,tag)
    coords= coords.reshape((-1,3))[:,0:2]


CreatingMesh(1,1,2,0.3,5)
GetData(1,1,2,0.3,5)
        