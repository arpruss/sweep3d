from struct import pack
import sys
import math
import cmath
import math
from vector import *
from triangulate import *

def toColorSCAD(polys, moduleName="object1"):
    def polyToSCAD(poly):
        pointsDict = {}
        i = 0
        out = "polyhedron(points=["
        points = []
        for face in poly:
            for v in face:
                if tuple(v) not in pointsDict:
                    pointsDict[tuple(v)] = i
                    points.append( "[%.9f,%.9f,%.9f]" % tuple(v) )
                    i += 1
        out += ",".join(points)
        out += "], faces=["
        out += ",".join( "[" + ",".join(str(pointsDict[tuple(v)]) for v in face) + "]" for face in poly ) + "]"
        out += ");"
        return out

    out = "module " +moduleName+ "() {\n"
    for rgb,monoPolys in polys:
        for poly in monoPolys:
            out += "  color([%.3f,%.3f,%.3f]) " % ( rgb[0]/255., rgb[1]/255., rgb[2]/255. ) 
            out += "%s\n" % polyToSCAD(poly)
    out += "}\n"
    return out
        
def saveSCAD(filename, polys, moduleName="object1", quiet=False):
    if not quiet: sys.stderr.write("Saving %s\n" % filename)
    with open(filename, "w") as f:
        f.write(toColorSCAD(polys, moduleName=moduleName))
        f.write("\n" + moduleName + "();")

def saveSTL(filename, mesh, swapYZ=False, quiet=False):
    if not quiet: sys.stderr.write("Saving %s\n" % filename)
    minY = float("inf")
    minVector = Vector(float("inf"),float("inf"),float("inf"))
    numTriangles = 0
    if swapYZ:
        matrix = Matrix( (1,0,0), (0,0,-1), (0,1,0) )
    else:
        matrix = Matrix.identity(3)
    for rgb,monoMesh in mesh:
        for triangle in monoMesh:
            numTriangles += 1
            for vertex in triangle:
                vertex = matrix*vertex
                minVector = Vector(min(minVector[i], vertex[i]) for i in range(3))
    minVector -= Vector(0.001,0.001,0.001) # make sure all STL coordinates are strictly positive as per Wikipedia
     
    with open(filename, "wb") as f:
        f.write(pack("80s",b''))
        f.write(pack("<I",numTriangles))
        for rgb,monoMesh in mesh:
            color = 0x8000 | ( (rgb[0] >> 3) << 10 ) | ( (rgb[1] >> 3) << 5 ) | ( (rgb[2] >> 3) << 0 )
            for tri in monoMesh:
                normal = (Vector(tri[1])-Vector(tri[0])).cross(Vector(tri[2])-Vector(tri[0])).normalize()
                f.write(pack("<3f", *(matrix*normal)))
                for vertex in tri:
                    f.write(pack("<3f", *(matrix*(vertex-minVector))))
                f.write(pack("<H", color))            
                
class SectionAligner(object):
    def __init__(self, upright=Vector(0,0,1), keepSectionUpright=False):
        self.upright = upright.normalize()
        self.keepSectionUpright = keepSectionUpright
        
        # perp will be perpendicular to upright, and serve as a default direction for the cross-section's normal
        m = min(abs(comp) for comp in upright)
        if self.upright.x == m:
            self.perp1 = (1,0,0) - self.upright.x*self.upright
        elif upright.y == m:
            self.perp1 = (0,1,0) - self.upright.y*self.upright
        else:
            self.perp1 = (0,0,1) - self.upright.z*self.upright
        self.perp2 = self.upright.cross(self.perp1)
        
    # First, the cross-section will be stood up with x-axis going to -perp2, y-axis going to upright and normal going to perp1
        self.m1 = Matrix( (-self.perp2.x, self.upright.x, 0), (-self.perp2.y, self.upright.y, 0), (-self.perp2.z, self.upright.z, 0) )
    
    def align(self, sectionPoints, direction, position):
        out = []
        direction = direction.normalize()
        projDirection = (direction - (direction*self.upright)*self.upright).normalize()
        
        # Then it will be rotated to match horizontal angle of the knot direction (which had better not be straight up or down along upright) 
        m2 = Matrix.rotateVectorToVector(self.perp1, projDirection)
        
        if self.keepSectionUpright:
            m = m2 * self.m1
        else:
            # Finally, we will tilt it up or down to match the direction
            m3 = Matrix.rotateVectorToVector(projDirection, direction)
            m = m3 * m2 * self.m1
        
        for v in sectionPoints:
            out.append( m * Vector(v) + position )

        return out
                
def sweep(mainPath, section, t1, t2, tstep, upright=Vector(0,0,1), 
        solid=False, clockwise=False, scad=False, cacheTriangulation=False, closed=True, keepSectionUpright=False):
    """
    The upright vector specifies the preferred pointing direction for the y-axis in the input sections.
    The tangent to the mainPain should never be close to parallel to the upright vector. E.g., for a mainly
    horizontal knot, the default (0,0,1) setting should work.
    
    cacheTriangulation optimizes in case all the section triangulations are the same. The typical case is where
    either the section is constant, or the different sections are affine transforms of one another.
    
    In polyhedron mode, each little piece of the knot is a separate polyhedron, appropriate for dumping into OpenSCAD.
    """
    
    if scad:
        solid = True
        clockwise = True
        
    upright = Vector(upright)
    cachedTriangulation = None

    aligner = SectionAligner(upright=upright, keepSectionUpright=keepSectionUpright)
    
    def getCrossSection(s, t):
        f1 = mainPath(t)
        if closed or t+tstep/2. <= t2:
            direction = mainPath( (t-t1+tstep/2.)%(t2-t1) + t1 ) - f1
        else:
            direction = f1 - mainPath( t-tstep/2. )
        return aligner.align(s, direction, f1)
        
    if (solid or not closed) and cacheTriangulation:
        cachedTriangulation = triangulate(section(t1))
        
    def getTriangulation(s):
        if cachedTriangulation is None:
            return triangulate(s)
        else:
            return cachedTriangulation
        
    output = []
    nextCrossSection = None
    t = t1
    while t < t2:
        if nextCrossSection is not None:
            curCrossSection,curTriangulation = nextCrossSection,nextTriangulation
        else:
            s1 = section(t)
            curCrossSection,curTriangulation = getCrossSection(s1, t), not solid or getTriangulation(s1)
        nextT = t+tstep if t+tstep < t2 else (t1 if closed else t2)
        s2 = section(nextT)
        nextCrossSection,nextTriangulation = getCrossSection( s2, nextT ), not solid or getTriangulation(s2)
        
        n = len(curCrossSection)
        assert n == len(nextCrossSection)
        
        def applyTriangulation(points,tr,clockwise=False):
            if clockwise:
                return [ (points[t[2]],points[t[1]],points[t[0]]) for t in tr ]
            else:
                return [ (points[t[0]],points[t[1]],points[t[2]]) for t in tr ]

        if not closed and not solid and t == t1:
            curTriangulation = getTriangulation(s1)
            output += applyTriangulation(curCrossSection, curTriangulation, clockwise=clockwise)
        
        def triangulateTube(i):
            if clockwise:
                return [ ( nextCrossSection[(i+1)%n], nextCrossSection[i], curCrossSection[i] ), 
                         ( curCrossSection[i], curCrossSection[(i+1)%n], nextCrossSection[(i+1)%n] ) ]
            else:
                return [ ( curCrossSection[i], nextCrossSection[i], nextCrossSection[(i+1)%n] ), 
                         ( nextCrossSection[(i+1)%n], curCrossSection[(i+1)%n], curCrossSection[i] ) ]
                     
        if solid:
            polyhedron = applyTriangulation(curCrossSection, curTriangulation, clockwise=clockwise)
            polyhedron += applyTriangulation(nextCrossSection, nextTriangulation, clockwise=not clockwise)
            for i in range(len(curCrossSection)):
                polyhedron += triangulateTube(i)
            output.append(polyhedron)
        else:
            for i in range(n):
                output += triangulateTube(i)
                
        t += tstep
            
    if not solid and not closed:
        s = section(t2)
        output += applyTriangulation(getCrossSection(s, t2), getTriangulation(s), clockwise=not clockwise)
        
    return output
    
if __name__ == '__main__':    
    r = math.sqrt(3)/3.
    scale = 5
    path1 = lambda t: scale*Vector( math.cos(t), math.sin(t)+r, -math.cos(3*t)/3.  )
    path2 = lambda t: scale*Vector( math.cos(t)+0.5, math.sin(t)-r/2., -math.cos(3*t)/3. )
    path3 = lambda t: scale*Vector( math.cos(t)-0.5, math.sin(t)-r/2, -math.cos(3*t)/3. )

    # we define the base section with complex numbers to make rotation simple
    #baseSection = Vector( cmath.exp(1j*2*math.pi*k/3) for k in range(3) ) #triangle
    #baseSection = Vector( 0j, 0.5+0j,0.5+.2j,0.2+.2j,0.2+1j,1j ) # letter L
    #baseSection = Vector( 1+0j,1j,-1+0j,0-1j ) #square
    baseSection = Vector( cmath.exp(2j*math.pi*k/10) * (1 if k%2 else 0.5) for k in range(10) ) #star with five points

    # now spin the baseSection as we sweep it
    section = lambda t : cmath.exp(1j*t) * baseSection

    rings = []
    rings.append( ( (255,0,0), sweep(path1, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True) ) )
    rings.append( ( (0,255,0), sweep(path2, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True) ) )
    rings.append( ( (0,0,255), sweep(path3, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True) ) )
#    rings.append( ( (255,0,0), sweep(lambda t:Vector(t,0,0), section, 0,6,.1,scad=True, closed=False)))

    saveSCAD("borromean.scad", rings)

    rings = []
    rings.append( ( (255,0,0), sweep(path1, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=False) ) )
    rings.append( ( (0,255,0), sweep(path2, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=False) ) )
    rings.append( ( (0,0,255), sweep(path3, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=False) ) )
#    rings.append( ( (255,0,0), sweep(lambda t:Vector(t,0,0), section, 0,6,.1,closed=False)))

    saveSTL("borromean.stl", rings)

    rings = []
    # cinquefoil from http://www.maa.org/sites/default/files/images/upload_library/23/stemkoski/knots/page6.html
    
    cinqueFoilPath = lambda t: scale/2.*Vector(math.cos(2*t) * (3 + math.cos(5*t)), math.sin(2*t) * (3 + math.cos(5*t)), math.sin(5*t))
    saveSTL("cinquefoil.stl", [ ( (255,255,0), sweep(cinqueFoilPath, section, 0, 2*math.pi, .05) ) ] )
    saveSCAD("cinquefoil.scad", [ ( (255,255,0), sweep(cinqueFoilPath, section, 0, 2*math.pi, .05, scad=True, cacheTriangulation=True) ) ] )
             
    # screw thread parametrized by number of turns
    screwLength = 25 # this is nominal: the shaft will stick out upwards to accommodate the last bit of thread
    shaftDiameter = 10
    pitch = 5
    nTurns = screwLength / float(pitch)
    taperLength = 0.4 # in turns
    threadBase = Vector( Vector(0,-0.5), Vector(0.5,0), Vector(0,0.5) )
    threadHeightPerPitch = 0.75 # fraction of pitch taken up by thread
    
    # draw thread
    def threadSection(t): # if you change it so it's not affine, you'll have to disable triangulation caching
        if t < taperLength:
            fraction = t/taperLength
            if fraction < 0.01:
                fraction = 0.01 # we don't want a degenerate section
        else:
            fraction = 1.
        return threadHeightPerPitch * pitch * fraction * threadBase
    def threadPath(t):
        angle = 2 * math.pi * t
        return Vector( 0.5 * shaftDiameter * math.cos(angle), 0.5 * shaftDiameter * math.sin(angle), t * screwLength / nTurns )
    screw = []
    screw.append( ( (0,0,255), sweep(threadPath, threadSection, 0, nTurns, 1./16, upright=Vector(0,0,1), 
            keepSectionUpright=True, closed=False, cacheTriangulation=True, scad=True)) )
    
    # this could be just a cylinder in OpenSCAD, but it's fun to show how to do it using sweep
    circularPrecision = 32
    adjDiameter = 1.0000001 * shaftDiameter / math.cos(math.pi/circularPrecision)
    adjHeight = screwLength + pitch * threadHeightPerPitch / 2
    baseSection = adjDiameter/2 * Vector( cmath.exp(2j * math.pi * k / circularPrecision) for k in range(circularPrecision) )
    shaftPath = lambda t : Vector(0, 0, t * adjHeight )
    
    screw.append( ( (0,0,128), sweep(shaftPath, lambda t:baseSection, 0, 1, 1, upright=(1,0,0), scad=True, keepSectionUpright=True, 
            closed=False, cacheTriangulation=True )))
    
    saveSCAD("screwthread.scad", screw)
        
