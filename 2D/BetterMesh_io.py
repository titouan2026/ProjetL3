import gmsh_api.gmsh as gmsh
import numpy as np
import meshio

def CreateMesh(R,nm,nd,lambda0,n,geo):
    
    lc=1
    c=np.real(nd)/nm
    
    global object
    global objectLines
    
    object=[]
    objectLines=[]
    
    gmsh.finalize()
    gmsh.initialize()
    
    class Object:
        def __init__(self, x1, y1, x2, y2, x3, y3, x4, y4):
            global object
            global objectLines
            self.x1 = x1
            self.y1 = y1
            self.x2 = x2
            self.y2 = y2
            self.x3 = x3
            self.y3 = y3
            self.x4 = x4
            self.y4 = y4
            
            self.point1 = gmsh.model.occ.add_point(self.x1, self.y1, 0, lc/n)
            self.point2 = gmsh.model.occ.add_point(self.x2, self.y2, 0, lc/n)
            self.point3 = gmsh.model.occ.add_point(self.x3, self.y3, 0, lc/n)
            self.point4 = gmsh.model.occ.add_point(self.x4, self.y4, 0, lc/n)

            self.line1 = gmsh.model.occ.add_line(self.point1, self.point2)
            self.line2 = gmsh.model.occ.add_line(self.point2, self.point3)
            self.line3 = gmsh.model.occ.add_line(self.point3, self.point4)
            self.line4 = gmsh.model.occ.add_line(self.point4, self.point1)
            
            self.curveloop = gmsh.model.occ.addCurveLoop([self.line1,self.line2,self.line3,self.line4])
            
            object.append(self.curveloop)
            objectLines+=[self.line1,self.line2,self.line3,self.line4]
            
            
            

    gmsh.option.setNumber("General.Verbosity", 2)
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 2*lambda0/n)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5*lambda0/(c*n))

    circle=gmsh.model.occ.addCircle(0.0, 0.0, 0.0,R,5,angle1=0.,angle2=2*np.pi)
    
    #création des objet en fonction de la liste de géométrie
    for obj in geo[0]:
        x1,y1,x2,y2,x3,y3,x4,y4=obj
        Object(x1, y1, x2, y2, x3, y3, x4, y4)

    circle_loop=gmsh.model.occ.addCurveLoop([circle])

    gmsh.model.occ.addPlaneSurface([circle_loop]+object,8)
    
    
    object_surfaces=[]
    
    for i, obj in enumerate(object):
        object_surfaces.append(gmsh.model.occ.addPlaneSurface([obj]))


    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(1,[circle],tag=16)
    gmsh.model.addPhysicalGroup(1,objectLines,tag=19)
    gmsh.model.addPhysicalGroup(2,object_surfaces,tag=17)
    gmsh.model.addPhysicalGroup(2,[8],tag=18)
    
    field_size = gmsh.model.mesh.field.add("MathEval",13)
    gmsh.model.mesh.field.setString(field_size, "F", "{}".format(lambda0/n))

    field_size2 = gmsh.model.mesh.field.add("MathEval",14)
    gmsh.model.mesh.field.setString(field_size2, "F", "{}".format(lambda0/(c*n)))

    field_restrict = gmsh.model.mesh.field.add("Restrict",15)
    gmsh.model.mesh.field.setNumbers(field_restrict, "SurfacesList", object_surfaces)
    gmsh.model.mesh.field.setNumber(field_restrict, "IField", field_size2)

    gmsh.model.mesh.field.add("Min", 12)
    gmsh.model.mesh.field.setNumbers(12, "FieldsList",[15,13])
    gmsh.model.mesh.field.setAsBackgroundMesh(12)

    gmsh.model.mesh.generate(2)

    gmsh.write("C:/Users/PC/Desktop/cours/ProjetL3/Meshes/{}_{}_{}_{}_{}_{}.msh".format(R, nm, nd, lambda0, n, geo[1]))
    
    ### DEBUGGING TOOL ###
    
    gmsh.fltk.run()
    
    ######################
    
    gmsh.finalize()

def GetData(R,nm,nd,lambda0,n,geo_name):
    
    dim=-1
    tag=-1
    gmsh.initialize()
    
    # gmsh.option.setNumber("General.Verbosity", 2)
    mesh = meshio.read("C:/Users/PC/Desktop/cours/ProjetL3/Meshes/{}_{}_{}_{}_{}_{}.msh".format(R, nm, nd, lambda0, n, geo_name))
    
    gmsh.finalize()
    return mesh