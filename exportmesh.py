from struct import pack
from vector import *
import sys

def toSCADModule(polys, moduleName):
    """
    INPUT:
    polys: list of (color,polyhedra) pairs (counterclockwise triangles)
    moduleName: OpenSCAD module name
    
    OUTPUT: string with OpenSCAD code implementing the polys
    """
    scad = []
    scad.append("module " +moduleName+ "() {")
    for rgb,poly in polys:
        if rgb is not None:
            line = "  color([%.3f,%.3f,%.3f]) " % tuple(min(max(c,0.),1.) for c in rgb)
        else:
            line = ""
        pointsDict = {}
        i = 0
        line += "polyhedron(points=["
        points = []
        for face in poly:
            for v in reversed(face):
                if tuple(v) not in pointsDict:
                    pointsDict[tuple(v)] = i
                    points.append( "[%.9f,%.9f,%.9f]" % tuple(v) )
                    i += 1
        line += ",".join(points)
        line += "], faces=["
        line += ",".join( "[" + ",".join(str(pointsDict[tuple(v)]) for v in reversed(face)) + "]" for face in poly ) + "]"
        line += ");"
        scad.append(line)
    scad.append("}")
    return "\n".join(scad)

def saveSCAD(filename, polys, moduleName="object1", quiet=False):
    """
    filename: filename to write OpenSCAD file
    polys: list of (color,polyhedra) pairs (counterclockwise triangles)
    moduleName: OpenSCAD module name
    quiet: give no status message if set
    """
    if not quiet: sys.stderr.write("Saving %s\n" % filename)
    with open(filename, "w") as f:
        f.write(toSCADModule(polys, moduleName))
        f.write("\n" + moduleName + "();\n")

def saveSTL(filename, mesh, swapYZ=False, quiet=False):
    """
    filename: filename to save STL file
    mesh: list of (color,triangle) pairs (counterclockwise)
    swapYZ: should Y/Z axes be swapped?
    quiet: give no status message if set
    """
    if not quiet: sys.stderr.write("Saving %s\n" % filename)
    minY = float("inf")
    minVector = Vector(float("inf"),float("inf"),float("inf"))
    numTriangles = 0
    if swapYZ:
        matrix = Matrix( (1,0,0), (0,0,-1), (0,1,0) )
    else:
        matrix = Matrix.identity(3)
        
    mono = True
    for rgb,triangle in mesh:
        if rgb is not None:
            mono = False
        numTriangles += 1
        for vertex in triangle:
            vertex = matrix*vertex
            minVector = Vector(min(minVector[i], vertex[i]) for i in range(3))
    minVector -= Vector(0.001,0.001,0.001) # make sure all STL coordinates are strictly positive as per Wikipedia
     
    with open(filename, "wb") as f:
        f.write(pack("80s",b''))
        f.write(pack("<I",numTriangles))
        for rgb,tri in mesh:
            if mono:
                color = 0
            else:
                if rgb is None:
                    rgb = (255,255,255)
                else:
                    rgb = tuple(min(255,max(0,int(0.5 + 255 * comp))) for comp in rgb)
                color = 0x8000 | ( (rgb[0] >> 3) << 10 ) | ( (rgb[1] >> 3) << 5 ) | ( (rgb[2] >> 3) << 0 )
            normal = (Vector(tri[1])-Vector(tri[0])).cross(Vector(tri[2])-Vector(tri[0])).normalize()
            f.write(pack("<3f", *(matrix*normal)))
            for vertex in tri:
                f.write(pack("<3f", *(matrix*(vertex-minVector))))
            f.write(pack("<H", color))            
                