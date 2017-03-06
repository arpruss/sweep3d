from exportmesh import *
from vector import *
from sys import argv
import random
import math

nIterations = 2 if len(argv) < 2 else int(argv[1])

height = lambda n : 10. / 2**n

def r(n):
    h = height(n)
    return random.uniform(-0.5*h,0.5*h)
    
center = r(0)
corners = [r(0) for i in range(6)]
theta = math.pi / 3.
x = [10.*math.cos(theta*k) for k in range(6)]
y = [10.*math.sin(theta*k) for k in range(6)]

mesh = [ (Vector(0,0,center), Vector(x[k],y[k],corners[k]), 
               Vector(x[(k+1)%6],y[(k+1)%6],corners[(k+1)%6])) for k in range(6)]

for iteration in range(1,1+nIterations):
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

offset = Vector(0,0,0.01-min(v.z for face in mesh for v in face))

def project(v):
    return Vector(v[0],v[1],0)
            
bottom = [ tuple(project(v) for v in reversed(face)) for face in mesh ]
top = [ tuple(v+offset for v in face) for face in mesh ]
edges = set()

for face in top:
    for i in range(3):
        j = (i+1)%3
        revEdge = (face[j],face[i])
        if revEdge in edges:
            edges.remove(revEdge)
        else:
            edges.add((face[i],face[j]))

for edge in edges:
    top.append((project(edge[0]),project(edge[1]),edge[0]))
    top.append((project(edge[1]),edge[1],edge[0]))

mesh = [ (None, face) for face in bottom+top ]

saveSTL("noise.stl", mesh)
