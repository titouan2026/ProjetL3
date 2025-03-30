import gmsh_api.gmsh as gmsh
import numpy as np
import meshio

def CreateMesh(L,l,nm,nd,lambda0,n):
    
    lc=1
    c=np.real(nd)/nm
    c= np.round(c,5)
    gmsh.finalize()
    gmsh.initialize()
    
    gmsh.option.setNumber("General.Verbosity", 2)

    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lambda0/n)

    
    point1 = gmsh.model.occ.add_point(0.6*L, l, 0, lc/c)
    point2 = gmsh.model.occ.add_point(0.6*L, -l, 0, lc/c)
    point3 = gmsh.model.occ.add_point(0.2*L, -l, 0, lc/c)
    point4 = gmsh.model.occ.add_point(0.2*L, l, 0, lc/c)
    
    point5 = gmsh.model.occ.add_point(L, l, 0, lc)
    point6 = gmsh.model.occ.add_point(L, -l, 0, lc)
    point7 = gmsh.model.occ.add_point(0, -l, 0, lc)
    point8 = gmsh.model.occ.add_point(0, l, 0, lc)

    line1 = gmsh.model.occ.add_line(point1, point2)
    line2 = gmsh.model.occ.add_line(point2, point3)
    line3 = gmsh.model.occ.add_line(point3, point4)
    line4 = gmsh.model.occ.add_line(point4, point1)
    
    line5 = gmsh.model.occ.add_line(point5, point6)
    line6 = gmsh.model.occ.add_line(point6, point2)
    line7 = gmsh.model.occ.add_line(point3, point7)
    line8 = gmsh.model.occ.add_line(point7, point8)
    line9 = gmsh.model.occ.add_line(point8, point4)
    line10 = gmsh.model.occ.add_line(point1, point5)

    gmsh.model.occ.synchronize()
    

    gmsh.model.occ.addCurveLoop([line10,line5,line6,line1] ,30)
    gmsh.model.occ.addCurveLoop([line8,line9,line3,line7] ,31)
    gmsh.model.occ.addCurveLoop([line5,line6,line2,line7,line8,line9,line4,line10] ,6)
    gmsh.model.occ.addCurveLoop([line1,line2,line3,line4] ,7)
    gmsh.model.occ.addPlaneSurface([30],8)
    gmsh.model.occ.addPlaneSurface([31],32)

    gmsh.model.occ.addPlaneSurface([7],9) #Suface du millieu, peut etre modifié mais dois resté numéo 9!

    gmsh.model.occ.synchronize()
    
    gmsh.model.addPhysicalGroup(1,[line5,line8],tag=16)
    gmsh.model.addPhysicalGroup(1,[line1,line2,line3,line4,line6,line7,line9,line10],tag=19)
    gmsh.model.addPhysicalGroup(2,[9],tag=17)
    gmsh.model.addPhysicalGroup(2,[8,32],tag=18)
    
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

    gmsh.write("C:/Users/PC/Desktop/cours/ProjetL3/Meshes/{}_{}_{}_{}_{}_{}.msh".format(L,l,nm,nd,lambda0,n))
    # gmsh.fltk.run()
    
    gmsh.finalize()

def GetData(L,l,nm,nd,lambda0,n):
    
    dim=-1
    tag=-1
    gmsh.initialize()
    
    # gmsh.option.setNumber("General.Verbosity", 2)
    mesh = meshio.read("C:/Users/PC/Desktop/cours/ProjetL3/Meshes/{}_{}_{}_{}_{}_{}.msh".format(L,l,nm,nd,lambda0,n))
    
    gmsh.finalize()
    return mesh
