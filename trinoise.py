from __future__ import division
from exportmesh import *
from vector import *
import sys
import random
import math
import os
import getopt

def help():
    print "python trinoise [options] [outname]"

nIterations = 2
size = 10
amplitude = None
mode = "h"
 
try:
    opts, args = getopt.getopt(sys.argv[1:], "h", 
                    ["amplitude=", "size=", "iterations=", "sphere", "triangle", "hex", "help"])

    i = 0
    while i < len(opts):
        opt,arg = opts[i]
        if opt in ('-h', '--help'):
            help()
            sys.exit(0)
        elif opt == "--iterations":
            nIterations = int(arg)
        elif opt == "--amplitude":
            amplitude = float(arg)
        elif opt == "--size":
            size = float(arg)
        elif opt == "--sphere":
            mode = "s"
        elif opt == "--triangle":
            mode = "t"
        elif opt == "--hexagon":
            mode = "h"
            
        i += 1
except getopt.GetoptError as e:
    sys.stderr.write(str(e)+"\n")
    help()
    sys.exit(2)
        
if len(args):
    outName = args[0]
    if not os.path.splitext(outName)[1]:
        outName += ".stl"
else:
    outName = "noise.stl"
        
if amplitude == None:
    amplitude = 1 if mode == "s" else 10
                
amp = lambda n : 1. / 2**(n+1)

def r(n):
    a = amp(n)
    return random.uniform(-0.5*a,0.5*a)
    
if mode == "h":
    theta = math.pi / 3.
    x = [math.cos(theta*k) for k in range(6)]
    y = [math.sin(theta*k) for k in range(6)]

    mesh = [ (Vector(0,0), Vector(x[k],y[k]), Vector(x[(k+1)%6],y[(k+1)%6])) for k in range(6)]
elif mode == "t":
    mesh = [ (Vector(0,0), Vector(1,0), Vector(0.5,math.sqrt(3)/2.)) ]
elif mode == "s":
    vv = (Vector(-1,0,-1./math.sqrt(2)), Vector(1,0,-1./math.sqrt(2)), Vector(0,-1,1./math.sqrt(2)), Vector(0,1,1./math.sqrt(2)))
    vv = tuple(v.normalize() for v in vv)
    mesh = [(vv[2],vv[0],vv[1]), (vv[1],vv[0],vv[3]), (vv[0],vv[2],vv[3]), (vv[2],vv[1],vv[3])]
        
noise = {}
vertices = set( v for face in mesh for v in face )

for v in vertices:
    if v not in noise:
        noise[v] = r(0)

for iteration in range(1,1+nIterations):
    newMesh = []
    refinements = {}
    
    def refine(a,b):
        index = (min(a,b),max(a,b))
        try:
            return refinements[index]
        except KeyError:
            v = 0.5 * (a+b)
            if mode == "s":
                v = v.normalize()
            refinements[index] = v
            noise[v] = r(iteration)+0.5*(noise[a]+noise[b])
            vertices.add(v)
            return v
    
    for face in mesh:
        v1,v2,v3 = face
        newMesh.append( (v1,refine(v1,v2),refine(v1,v3))) 
        newMesh.append( (refine(v1,v2), v2, refine(v2,v3)))
        newMesh.append( (refine(v2,v3), v3, refine(v1,v3)))
        newMesh.append( (refine(v2,v3), refine(v1,v3), refine(v1,v2)))
    
    mesh = newMesh

if mode == "s":
    outMesh = []
    def adjust(v):
        newRadius = size + amplitude * noise[v]
        if newRadius < 0.01:
            newRadius = 0.01
        return newRadius * v
    for face in mesh:
        outMesh.append( (None, tuple( adjust(v) for v in face )) )
else:
    outMesh = []
    lowest = min(noise[v] for v in vertices)
    for v in vertices:
        noise[v] = amplitude * (noise[v]-lowest+0.01)
    bottom = [ tuple(size*v+Vector(0.,0.,0.) for v in reversed(face)) for face in mesh ]
    top = [ tuple(size*v+Vector(0.,0.,noise[v]) for v in face) for face in mesh ]
    
    edges = set()

    for face in top:
        for i in range(3):
            j = (i+1)%3
            revEdge = (face[j],face[i])
            if revEdge in edges:
                edges.remove(revEdge)
            else:
                edges.add((face[i],face[j]))
                
    def project(a):
        return Vector(a.x, a.y, 0.)

    for edge in edges:
        top.append((project(edge[0]),project(edge[1]),edge[0]))
        top.append((project(edge[1]),edge[1],edge[0]))
        
    outMesh = [ (None, face) for face in bottom+top ]

saveSTL(outName, outMesh)
