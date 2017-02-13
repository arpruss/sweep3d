from vector import *
from exportmesh import *
import itertools

def surfaceToMesh(data, center=False, twoSided=False, xScale=1., yScale=1., zScale=1., tolerance=1e-6, color=None):
    width = len(data)
    height = len(data[0])
    xMin = width - 1
    xMax = 0
    yMin = height - 1
    yMax = 0
        
    def getValue(*args):
        if len(args) == 1:
            x,y = args[0]
        else:
            x,y = args
            
        if x < 0 or x >= width or y < 0 or y >= height:
            return 0.
        else:
            if data[x][y] <= tolerance:
                return 0.
            else:
                return data[x][y]

    for x in range(width):
        for y in range(height):
            if getValue(x,y) > 0.: 
                xMin = min(xMin, x)
                xMax = max(xMax, x)
                yMin = min(yMin, y)
                yMax = max(yMax, y)
    
    if center:
        offsetX = -0.5 * (xMin + xMax)
        offsetY = -0.5 * (yMin + yMax)
    else:
        offsetX = - ( xMin - 1. )
        offsetY = - ( yMin - 1. )
        
    mesh = []
    
    for x in range(xMin - 1, xMax + 1):
        for y in range(yMin - 1, yMax + 1):
            v = Vector(x,y)
            numPoints = sum(1 for delta in ((0,0), (1,0), (0,1), (1,1)) if getValue(v+delta) > 0.)

            def triangles(d1, d2, d3):
                v1,v2,v3 = v+d1,v+d2,v+d3
                z1,z2,z3 = map(getValue, (v1,v2,v3))
                if (z1,z2,z3) == (0.,0.,0.):
                    return []
                z1,z2,z3 = zScale * Vector(z1,z2,z3)
                v1,v2,v3 = map((lambda w : Vector((w.x+offsetX)*xScale, (w.y+offsetY)*yScale)), (v1,v2,v3))
                output = [(color,((v1.x,v1.y,z1), (v2.x,v2.y,z2), (v3.x,v3.y,z3)))]
                if not twoSided:
                    z1,z2,z3 = 0.,0.,0.
                output.append ( (color,((v3.x,v3.y,-z3), (v2.x,v2.y,-z2), (v1.x,v1.y,-z1))) )
                return output
            
            if numPoints > 0:
                if getValue(v+(0,0)) == 0. and getValue(v+(1,1)) == 0.:
                    mesh += triangles((0,0), (1,0), (1,1))
                    mesh += triangles((1,1), (0,1), (0,0))
                else:
                    mesh += triangles((0,0), (1,0), (0,1))
                    mesh += triangles((1,0), (1,1), (0,1))

    return mesh
    
if __name__ == '__main__':
    data = [ [0.1,0.2], [0.3,0.4] ]
    saveSTL("surface.stl", surfaceToMesh(data, twoSided=True))