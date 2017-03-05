from exportmesh import *
from vector import *
from sys import argv

def moreFaces(v1,v2,v3):
    m = (1./3)*(v1+v2+v3)
    a = v2-v1
    b = v3-v1
    h = math.sqrt(6)/3. * a.norm()
    v4 = m + h * a.cross(b).normalize()
    return [ (v1,v2,v4), (v2,v3,v4), (v3,v1,v4) ]
    
v1,v2,v3 = Vector(0,0,0), Vector(10,math.sqrt(3)*10,0), Vector(20,0,0)

nIterations = 2 if len(argv) < 2 else int(argv[1])

mesh = [ (0,face) for face in moreFaces(v3,v2,v1)] + [(0,(v1,v2,v3))] 

for iteration in range(1,1+nIterations):
    newMesh = []
    for i,face in mesh:
        v1,v2,v3 = face
        newMesh.append((i,(v1,0.5*(v1+v2),0.5*(v1+v3))))
        newMesh.append((i,(0.5*(v1+v2), v2, 0.5*(v2+v3))))
        newMesh.append((i,(0.5*(v2+v3), v3, 0.5*(v1+v3))))
        newMesh += [(i+1,f) for f in moreFaces(0.5*(v2+v3), 0.5*(v1+v3), 0.5*(v1+v2))]
    mesh = newMesh
        
mesh = [ ((i/float(nIterations),0,1.-i/float(nIterations)), face) for i,face in mesh ]

saveSTL("fractal.stl", mesh)
