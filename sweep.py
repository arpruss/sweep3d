from struct import pack
import sys
import math
import cmath
import math
from vector import *
from triangulate import *

def toSCADModule(polys, moduleName):
    """
    INPUT:
    polys: list of (color,polyhedra) pairs (clockwise triangles)
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
            for v in face:
                if tuple(v) not in pointsDict:
                    pointsDict[tuple(v)] = i
                    points.append( "[%.9f,%.9f,%.9f]" % tuple(v) )
                    i += 1
        line += ",".join(points)
        line += "], faces=["
        line += ",".join( "[" + ",".join(str(pointsDict[tuple(v)]) for v in face) + "]" for face in poly ) + "]"
        line += ");"
        scad.append(line)
    scad.append("}")
    return "\n".join(scad)

def saveSCAD(filename, polys, moduleName="object1", quiet=False):
    """
    filename: filename to write OpenSCAD file
    polys: list of (color,polyhedra) pairs (clockwise triangles)
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
        solid=False, clockwise=False, scad=False, cacheTriangulation=False, closed=True, keepSectionUpright=False,
        color=None, tangentPrecision=0.5):
    """
    INPUTS:
    mainPath: function of t specifying the 3D path to sweep along -- output can be a vector or tuple
    section: either a list or a function of t returning a list, with contents being 2D points to sweep along the mainPath;
        the number of points must be the same for all values of t, and close values of t should result in similar
        sections
    t1,t2: start and end of t parameter for mainPath
    tstep: amount to increment t per step
    upright: vector to align the y-axis of the section towards; if the tangent of the mainPath ever aligns with upright,
        bad things will happen (numerical instability or code crash); ideally, you should point the upright vector away
        from the approximate plane of the mainPath curve
    solid: if True, returns a list of (color,polyhedron) pairs, with polyhedra composed of triangular faces; the polyhedra
        then should be joined together with a CSG program like OpenSCAD;
        if False, returns a list of (color,triangle) pairs for a mesh
    clockwise: if True, the order in a triangular face is clockwise when viewed from outside a solid
    scad: if True, sets solid=True and clockwise=True for OpenSCAD output
    cacheTriangulation: if True, it's assumed that all values of the section function can be triangulated in the same way; 
        the easiest way to guarantee this is to ensure that the different values of the section function are affine transforms
        of each other
    closed: specify whether the mainPath is a closed curve or not
    keepSectionUpright: if True, the y-axis of the section is exactly pointed along the upright vector; if False, the y-axis
        of the section is only pointed along the upright vector when the mainPath moves in the plane orthogonal to the upright
        vector; otherwise, the y-axis of the section gets tilted up or down to ensure that the section is normal to the
        mainPath's tangent
    color: either an RGB color (tuple/list of floats from 0 to 1) or a function from t to RGB colors (allowing a sweep varying
        in color); you can also specify None for the color, which if consistently done will generate a monochrome OpenSCAD file
        or STL mesh
    tangentPrecision: fraction of tstep used for calculating tangents
    
    OUTPUT:
    list of (color,polyhedron) or (color,triangle) pairs
    """
    
    if not hasattr(section, '__call__'):
        section0 = section
        section = lambda t : section0
        cacheTriangulation = True
    
    if not hasattr(color, '__call__'):
        color0 = color
        color = lambda t : color0
    
    if scad:
        solid = True
        clockwise = True
        
    upright = Vector(upright)
    cachedTriangulation = None

    aligner = SectionAligner(upright=upright, keepSectionUpright=keepSectionUpright)
    
    def getCrossSection(s, t):
        f1 = Vector(mainPath(t))
        if closed or t+tstep*tangentPrecision <= t2:
            direction = Vector(mainPath( (t-t1+tstep*tangentPrecision)%(t2-t1) + t1 )) - f1
        else:
            direction = f1 - Vector(mainPath( t-tstep*tangentPrecision ))
        return aligner.align(s, direction, f1)
        
    if (solid or not closed) and cacheTriangulation:
        cachedTriangulation = triangulate(section(t1))
        
    def getTriangulation(s):
        if cachedTriangulation is None:
            return triangulate(s)
        else:
            return cachedTriangulation
            
    curColor = None
    def colorize(data):
        return [(curColor, p) for p in data]
        
    output = []
    nextCrossSection = None
    t = t1
    while t < t2:
        curColor = color(t)
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
            output += colorize(applyTriangulation(curCrossSection, curTriangulation, clockwise=clockwise))
        
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
            output.append((curColor,polyhedron))
        else:
            for i in range(n):
                output += colorize(triangulateTube(i))
                
        t += tstep
            
    if not solid and not closed:
        s = section(t2)
        curColor = color(t2)
        output += colorize(applyTriangulation(getCrossSection(s, t2), getTriangulation(s), clockwise=not clockwise))
        
    return output
    
def scadScrew(screwLength, shaftDiameter, pitch, threadHeightPerPitch = 0.75, 
        threadBase = ( (0,-0.5), (0.5,0), (0,0.5) ), 
        upright=(0,0,1), start=(0,0,0), moduleName="screw", resolution=32, tolerance=0., leftHanded=False, clip=True):
    """
    INPUTS:
    screwLength: length of screw shaft
    shaftDiameter: diameter of screw shaft
    pitch: distance along screw between successive threads
    threadHeightPerPitch: scale of thread cross-section relative to pitch; if the threadBase has height 1,
        and threadHeigthPerPitch is 1, then the thread fills all of the pitch
    threadBase: thread profile, the default being a 45-90-45 triangle
    upright: direction of screw shaft
    start: position of start of screw shaft
    moduleName: name of OpenSCAD module containing this screw
    resolution: number of points per turn
    tolerance: thickens screw for subtracting from solids for female thread
    leftHanded: if True, thread is left-handed (threadBase is automatically rotated by 180 degrees)
    clip: if True, the thread doesn't stick out past the end of the shaft; of course, if making a real screw you generally want
        that; but if you are just generating a screw shape for subtraction from a solid for female through-thread, then things 
        will work faster without clipping
        
    OUTPUT:
    OpenSCAD module to generate the screw; if clipping is True, there will also be a moduleName_unclipped module generated
    """    
    
    nTurns = screwLength / float(pitch)
    
    if leftHanded:
        threadBase = Vector( -Vector(v) for v in threadBase )
        sign = -1.
    else:
        threadBase = Vector( Vector(v) for v in threadBase )
        sign = 1.
        
    upright = Vector(upright).normalize() 
    
    if tolerance != 0:
        threadHeightPerPitch = (threadHeightPerPitch * pitch + 2. * tolerance) / pitch
        toleranceInTurns = tolerance / float(pitch)
    else:
        toleranceInTurns = 0.
    
    def threadPath(t):
        angle = 2 * math.pi * t
        return Vector( 0.5 * shaftDiameter * math.cos(angle), 0.5 * sign * shaftDiameter * math.sin(angle), t * screwLength / nTurns ) + start

    screw = []
    screw += sweep(threadPath, threadHeightPerPitch * pitch * threadBase, -0.5-toleranceInTurns, nTurns+0.5+toleranceInTurns, 1./resolution, 
            upright=upright, keepSectionUpright=True, closed=False, cacheTriangulation=False, scad=True)
    
    # this could be just a cylinder in OpenSCAD, but it's fun to show how to do it using sweep
    adjDiameter = 1.001 * (shaftDiameter + 2 * tolerance) / math.cos(math.pi/resolution)
    baseSection = adjDiameter/2 * Vector( cmath.exp(2j * math.pi * k / resolution) for k in range(resolution) )

    shaftPath = lambda t : start - tolerance * upright + t * (screwLength + 2 * tolerance) * upright
    
    screw += sweep(shaftPath, baseSection, 0, 1, 1, upright=Vector(upright).perpendicular(), scad=True, 
            keepSectionUpright=True, closed=False, cacheTriangulation=True)
            
    screwSCAD = toSCADModule(screw, moduleName+("_unclipped" if clip else ""))
    
    if clip:
        screwSCAD += """

module %s() {
  render(convexity=1) difference() {
    %s_unclipped();
    translate([0,0,%.6f]) linear_extrude(height=%.6f) square(%.6f,center=true);
    translate([0,0,%.6f]) linear_extrude(height=%.6f) square(%.6f,center=true);
  }
}
""" % (moduleName, moduleName, -pitch,pitch,shaftDiameter+pitch*10, screwLength,pitch,shaftDiameter+pitch*10)

    return screwSCAD
    
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
    rings += sweep(path1, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True, color=(1.,0.,0.))
    rings += sweep(path2, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True, color=(0.,1.,0.))
    rings += sweep(path3, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True, color=(0.,0.,1.))

    saveSCAD("borromean.scad", rings)

    rings = []
    rings += sweep(path1, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=False, color=(1.,0.,0.))
    rings += sweep(path2, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=False, color=(0.,1.,0.))
    rings += sweep(path3, section, 0, 2*math.pi, .1, upright=Vector(0,.1,1), scad=False, color=(0.,0.,1.))

    saveSTL("borromean.stl", rings)

    rings = []
    # cinquefoil from http://www.maa.org/sites/default/files/images/upload_library/23/stemkoski/knots/page6.html
    
    cinqueFoilPath = lambda t: scale/2.*Vector(math.cos(2*t) * (3 + math.cos(5*t)), math.sin(2*t) * (3 + math.cos(5*t)), math.sin(5*t))
    color = lambda t: (1.,0.5*(math.cos(t)+1),0.)
    saveSTL("cinquefoil.stl", sweep(cinqueFoilPath, section, 0, 2*math.pi, .05, color=color))
    saveSCAD("cinquefoil.scad", sweep(cinqueFoilPath, section, 0, 2*math.pi, .05, scad=True, cacheTriangulation=True, color=color))
             
    screw = scadScrew(25, 10, 5, threadHeightPerPitch=0.75, resolution=40, moduleName="screw")

    screw += scadScrew(25, 10, 5, threadHeightPerPitch=0.75, resolution=40, moduleName="oversize_screw", tolerance=0.75, clip=False)
    screw += """

module nut() {
    render(convexity=1)
    difference() {
        cylinder(d=30,h=10,$fn=6);
        translate([0,0,0]) oversize_screw();
    }
}
    
translate([0,0,10-0.001]) screw();
cylinder(d=30,h=10,$fn=6);
translate([35,0,0]) nut();
"""
    
    sys.stderr.write("Saving screw.scad\n")
    with open("screw.scad", "wb") as f: f.write(screw)
    
