from exportmesh import *
from vector import *
from sys import argv
from random import uniform

nIterations = 2 if len(argv) < 2 else int(argv[1])

height = lambda n : 10. / 2**n

def r(n):
    h = height(n)
    return uniform(-0.5*h,0.5*h)

mesh = [ ( Vector(0,0,r(0)), Vector(10,math.sqrt(3)*10,r(0)), Vector(20,0,r(0)) ) ]

for iteration in range(nIterations):
    newMesh = []
    newHeights = {}
    
    def adjust(a,b):
        index = (min(a,b),max(a,b))
        try:
            return newHeights[index]
        except KeyError:
            v = 0.5 * (a+b) + (0,0,r(iteration))
            newHeights[index] = v
            return v
    
    for face in mesh:
        v1,v2,v3 = face
        newMesh.append( (v1,adjust(v1,v2),adjust(v1,v3))) 
        newMesh.append( (adjust(v1,v2), v2, adjust(v2,v3)))
        newMesh.append( (adjust(v2,v3), v3, adjust(v1,v3)))
        newMesh.append( (adjust(v2,v3), adjust(v1,v3), adjust(v1,v2)))
    
    mesh = newMesh
        
den = 1. if nIterations == 0 else float(nIterations)

mesh = [ (None, face) for face in mesh ]

saveSTL("noise.stl", mesh)
