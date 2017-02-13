from struct import pack
import sys
import cmath
import math
from vector import *
from triangulate import *
from exportmesh import *

class SectionAligner(object):
    def __init__(self, upright=Vector(0,0,1), keepSectionUpright=False):
        self.keepSectionUpright = keepSectionUpright
        self.upright = upright
        if not hasattr(self.upright, '__call__'):
            self.updateUpright(upright)
            self.uprightFunction = None
        else:
            self.uprightFunction = self.upright
    
    def updateUpright(self, upright):
        self.upright = Vector(upright).normalize()
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
    
    def align(self, sectionPoints, direction, position, t):
        out = []
        direction = direction.normalize()
        
        if self.uprightFunction is not None:
            self.updateUpright(self.uprightFunction(t))
        
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
        
def sweep(mainPath, section, t1, t2, steps=64, upright=Vector(0,0,1), derivative=None,
        solid=False, scad=False, cacheTriangulation=False, closed=True, keepSectionUpright=False,
        color=None, differentiationPrecision=0.05, symmetricDifferentiation=False):
    """
    INPUTS:
    mainPath: function of t specifying the 3D path to sweep along -- output can be a vector or tuple
    section: either a list or a function of t returning a list, with contents being 2D points to sweep along the mainPath;
        the number of points must be the same for all values of t, and close values of t should result in similar
        sections
    t1,t2: start and end of t parameter for mainPath
    steps: number of segments to include in sweep
    upright: upright vector or function of t returning an upright vector. The upright vector is what the y-axis of the
        section is aligned towards; if the tangent of the mainPath ever aligns with the upright vector,
        bad things will happen (numerical instability or code crash); ideally, you should point the upright vector away
        from the approximate plane of the mainPath curve; if you set the upright to None, the binormal will be calculated
        numerically and used.
    derivative: constant vector or function returning a vector. If None, this is calculated numerically. The derivative of the 
        curve is used as the normal for aligning the section. For special effects, you can set this to something other than
        the actual derivative (but note that this is what will be used to calculate the binormal if upright is None).
    solid: if True, returns a list of (color,polyhedron) pairs, with polyhedra composed of triangular faces; the polyhedra
        then should be joined together with a CSG program like OpenSCAD;
        if False, returns a list of (color,triangle) pairs for a mesh
    scad: if True, sets solid=True for OpenSCAD output; for backwards compatibility
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
    differentiationPrecision: will be scaled by tstep size
    symmetricDifferentiation: *TODO* currently produces bad screw thread
    
    OUTPUT:
    list of (color,polyhedron) or (color,triangle) pairs
    """
    
    steps = int(steps)
    tstep = (t2-t1)/float(steps)
    
    def ensureFunction(f):
        if hasattr(f, '__call__'):
            return f
        else: 
            return lambda t,v=f: v
                    
    if not hasattr(section, '__call__'):
        cachedTriangulation = True
        section = ensureFunction(section)
    color = ensureFunction(color)
    
    def differentiate(path, t, dt):
        if symmetricDifferentiation:
            if closed:
                def wrap(u):
                    return (u-t1) % (t2-t1) + t1
                sample1 = wrap(t-0.5*dt)
                sample2 = wrap(t+0.5*dt)
            elif t - 0.5 * dt < t1:
                sample1 = t1
                sample2 = t1 + dt
            elif t + 0.5 * dt > t2:
                sample1 = t2 - dt
                sample2 = t2
            else:
                sample1 = t - 0.5 * dt
                sample2 = t + 0.5 * dt
        else:
            if closed:
                def wrap(u):
                    return (u-t1) % (t2-t1) + t1
                sample1 = t
                sample2 = wrap(t+dt)
            elif t + dt > t2:
                sample1 = t2 - dt
                sample2 = t2
            else:
                sample1 = t
                sample2 = t + dt
        
        return Vector(path(sample2))-Vector(path(sample1))
    
    if derivative is None: 
        derivative = lambda t : differentiate(mainPath, t, 0.5*differentiationPrecision*tstep)
    else:
        derivative = ensureFunction(derivative)
        
    if upright is None:
        upright = lambda t: derivative(t).cross(differentiate(derivative,t,differentiationPrecision*tstep))
    
    if scad:
        solid = True
        
    cachedTriangulation = None

    aligner = SectionAligner(upright=upright, keepSectionUpright=keepSectionUpright)
    
    def getCrossSection(s, t):
        f1 = Vector(mainPath(t))
        return aligner.align(s, derivative(t), f1, t)
        
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
    
    for step in range(steps):
        t = t1 + step*(t2-t1)/float(steps)

        curColor = color(t)
        if nextCrossSection is not None:
            curCrossSection,curTriangulation = nextCrossSection,nextTriangulation
        else:
            s1 = section(t)
            curCrossSection,curTriangulation = getCrossSection(s1, t), not solid or getTriangulation(s1)

        if closed:
            nextT = (step+1)%steps * (t2-t1) / float(steps) + t1
        else:
            nextT = (step+1) * (t2-t1) / float(steps)

        s2 = section(nextT)
        nextCrossSection,nextTriangulation = getCrossSection( s2, nextT ), not solid or getTriangulation(s2)
        
        n = len(curCrossSection)
        assert n == len(nextCrossSection)
        
        def applyTriangulation(points,tr,reversed=False):
            if reversed:
                return [ (points[t[2]],points[t[1]],points[t[0]]) for t in tr ]
            else:
                return [ (points[t[0]],points[t[1]],points[t[2]]) for t in tr ]

        if not closed and not solid and t == t1:
            curTriangulation = getTriangulation(s1)
            output += colorize(applyTriangulation(curCrossSection, curTriangulation))
        
        def triangulateTube(i):
            return [ ( curCrossSection[i], nextCrossSection[i], nextCrossSection[(i+1)%n] ), 
                     ( nextCrossSection[(i+1)%n], curCrossSection[(i+1)%n], curCrossSection[i] ) ]
                     
        if solid:
            polyhedron = applyTriangulation(curCrossSection, curTriangulation)
            polyhedron += applyTriangulation(nextCrossSection, nextTriangulation, reversed=True)
            for i in range(len(curCrossSection)):
                polyhedron += triangulateTube(i)
            output.append((curColor,polyhedron))
        else:
            for i in range(n):
                output += colorize(triangulateTube(i))
            
    if not solid and not closed:
        s = section(t2)
        curColor = color(t2)
        output += colorize(applyTriangulation(getCrossSection(s, t2), getTriangulation(s)), reversed=True)
        
    return output
    
def regularPolygon(r=1,n=16,center=Vector(0.,0.)):
    return [r*Vector(math.cos(2*math.pi*k/n),math.sin(2*math.pi*k/n))+center for k in range(n)]
    
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
        
        TODO: clip right now only works for upright screws starting at (0,0,0) 
        
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
    screw += sweep(threadPath, threadHeightPerPitch * pitch * threadBase, -0.5-toleranceInTurns, nTurns+0.5+toleranceInTurns, steps=int(resolution*nTurns), 
            upright=upright, keepSectionUpright=True, closed=False, cacheTriangulation=False, scad=True)
    
    # this could be just a cylinder in OpenSCAD, but it's fun to show how to do it using sweep
    adjDiameter = 1.001 * (shaftDiameter + 2 * tolerance) / math.cos(math.pi/resolution)

    baseSection = regularPolygon(r=adjDiameter/2, n=resolution)

    shaftPath = lambda t : start - tolerance * upright + t * (screwLength + 2 * tolerance) * upright
    
    screw += sweep(shaftPath, baseSection, 0, 1, steps=1, upright=Vector(upright).perpendicular(), scad=True, 
            keepSectionUpright=True, closed=False, cacheTriangulation=True)
            
    screwSCAD = toSCADModule(screw, moduleName+("_unclipped" if clip else ""))
    
    if clip:
        screwSCAD += """

module %s() {
  render(convexity=1) intersection() {
    %s_unclipped();
    linear_extrude(height=%.6f) square(%.6f,center=true);
  }
}
""" % (moduleName, moduleName, screwLength,shaftDiameter+pitch*10)

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

    
    # Borromean rings
    rings = []
    rings += sweep(path1, section, 0, 2*math.pi, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True, color=(1.,0.,0.))
    rings += sweep(path2, section, 0, 2*math.pi, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True, color=(0.,1.,0.))
    rings += sweep(path3, section, 0, 2*math.pi, upright=Vector(0,.1,1), scad=True, cacheTriangulation=True, color=(0.,0.,1.))

    saveSCAD("borromean.scad", rings)

    rings = []
    rings += sweep(path1, section, 0, 2*math.pi, upright=Vector(0,.1,1), scad=False, color=(1.,0.,0.))
    rings += sweep(path2, section, 0, 2*math.pi, upright=Vector(0,.1,1), scad=False, color=(0.,1.,0.))
    rings += sweep(path3, section, 0, 2*math.pi, upright=Vector(0,.1,1), scad=False, color=(0.,0.,1.))

    saveSTL("borromean.stl", rings)

    # Cinquefoil from http://www.maa.org/sites/default/files/images/upload_library/23/stemkoski/knots/page6.html
    rings = []
    
    cinqueFoilPath = lambda t: scale/2.*Vector(math.cos(2*t) * (3 + math.cos(5*t)), math.sin(2*t) * (3 + math.cos(5*t)), math.sin(5*t))
    color = lambda t: (1.,0.5*(math.cos(t)+1),0.)
    saveSTL("cinquefoil.stl", sweep(cinqueFoilPath, section, 0, 2*math.pi, steps=128, color=color))
    saveSCAD("cinquefoil.scad", sweep(cinqueFoilPath, section(0), 0, 2*math.pi, steps=128, scad=True, cacheTriangulation=True, color=color))
    
    
       
    # Screw and nut (OpenSCAD only)
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
    with open("screw.scad", "w") as f: f.write(screw)
    
    # Threaded circle
    baseRadius = 40
    secondaryRadius = 10
    nThreads = 20
    threadSize = 0.75 * 2 * math.pi * baseRadius / nThreads
    
    def circularSpiral(t):
        r = secondaryRadius*math.cos(t*nThreads) + baseRadius
        z = -secondaryRadius*math.sin(t*nThreads)
        return Vector(r * math.cos(t), r * math.sin(t), z)

    threadSection = threadSize * Vector( Vector(0,-0.5), Vector(0.5,0), Vector(0,0.5) )
    
    # align the upright vector to be tangent to the base circle
    upright = lambda t:Vector(-math.sin(t),math.cos(t),0)
    
    color = lambda t:(0.5+0.5*math.cos(t),0.25,0.5+0.5*math.sin(t))
    
    threadedCircle = sweep(circularSpiral, threadSection, 0, 2*math.pi, steps=64*nThreads, scad=True, cacheTriangulation=True, closed=True,
                        upright=upright, keepSectionUpright=False, color=color)
                        
    threadedCircle += sweep(lambda t:(baseRadius*math.cos(t), baseRadius*math.sin(t), 0), regularPolygon(r=secondaryRadius*1.05), 0, 2*math.pi, 
                        scad=True, cacheTriangulation=True, closed=True, color=color)
    
    saveSCAD("threadedcircle.scad", threadedCircle)
