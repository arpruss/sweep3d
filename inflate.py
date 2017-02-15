from __future__ import division
from surface import *
import svgpath.shader as shader
import svgpath.parser as parser
import sys
import getopt
from exportmesh import *

def rasterizePolygon(polygon, spacing, shadeMode=shader.Shader.MODE_EVEN_ODD):
    """
    Returns boolean raster of strict interior as well as coordinates of lower-left corner.
    """
    spacing = float(spacing)
    lines = shader.Shader.shadePolygon(polygon, 0, spacing, avoidOutline=True, mode=shadeMode, alternate=False)
    height = len(lines)
    bottom = min(z[0].imag for z in lines)
    left = min(z[0].real for z in lines) - 0.5 * spacing
    right = max(z[1].real for z in lines)
    width = int((right-left) / spacing+1)

    raster = [[False for y in range(height)] for x in range(width)]

    for line in lines:
        y = int((line[0].imag - bottom) / spacing)
        x = int((line[0].real - left) / spacing)
        while left + x * spacing < line[1].real:
            if line[0].real < left + x * spacing:
                raster[x][y] = True
            x += 1
    
    return raster,complex(left,bottom)
    
def inflatePolygon(polygon, spacing=1., shadeMode=shader.Shader.MODE_EVEN_ODD, thickness=10., roundness=1., iterations=None, center=False, twoSided=False):
    # polygon is described by list of (start,stop) pairs, where start and stop are complex numbers
    raster,bottomLeft = rasterizePolygon(polygon, spacing, shadeMode=shadeMode)
    rasterWidth = len(raster)
    rasterHeight = len(raster[0])
    bottomLeftV = Vector(bottomLeft)
    surface = inflateRaster(raster, thickness=thickness, roundness=roundness, iterations=iterations)
    mesh0 = surfaceToMesh(surface, center=False, twoSided=twoSided)
    
    def scaleFace(face):
        return tuple(v*spacing+bottomLeftV for v in face)

    def trimFace(face, polygon, shadeMode=shader.Shader.MODE_EVEN_ODD):
        def trimLine(start, stop):
            delta = (stop - start).toComplex() # projects to 2D
            if delta == 0j:
                return stop
            length = abs(delta)
            normDelta = delta / length
            rotate = 1. / normDelta
            z0 = start.toComplex()
            
            class State(object): pass
            state = State()
            state.changed = False
            state.bestLength =  length
            
            for line in polygon:
                def update(x):
                    if 0 <= x < state.bestLength:
                        state.changed = True
                        state.bestLength = x
            
                l0 = rotate * (line[0]-z0)
                l1 = rotate * (line[1]-z0)
                if l0.imag == l1.imag and l0.imag == 0.:
                    if l0.real <= 0 and l1.real >= 0:
                        return start
                    update(l0.real)
                    update(l1.real)
                elif l0.imag <= 0 <= l1.imag or l1.imag <= 0 <= l0.imag:
                    # crosses real line
                    mInv = (l1.real-l0.real)/(l1.imag-l0.imag)
                    # (x - l0.real) / mInv = y - l0.imag
                    # so for y = 0: 
                    x = -l0.imag * mInv + l0.real
                    update(x)
            if state.changed:
                z = z0 + state.bestLength * normDelta
                return Vector(z.real, z.imag, 0.)
            else:
                print "not found"
                return stop
    
        def inside(v):
            if v[0] < 0 or v[0] >= rasterWidth or v[1] < 0 or v[1] >= rasterHeight:
                return False
            return raster[v[0]][v[1]]
            
        outsideCount = sum(1 for v in face if not inside(v))
        if outsideCount == 3 or outsideCount == 0:
            return [scaleFace(face)]
        elif outsideCount == 2:
            if inside(face[1]):
                face = (face[1], face[2], face[0])
            elif inside(face[2]):
                face = (face[2], face[0], face[1])
            # now, the first vertex is inside and the others are outside
            face = scaleFace(face)
            return [ (face[0], trimLine(face[0], face[1]), trimLine(face[0], face[2])) ]
        else: # outsideCount == 1
            if not inside(face[0]):
                face = (face[1], face[2], face[0])
            elif not inside(face[1]):
                face = (face[2], face[0], face[1])
            # now, the first two vertices are inside, and the third is outside
            face = scaleFace(face)
            closest0 = trimLine(face[0], face[2])
            closest1 = trimLine(face[1], face[2])
            if closest0 != closest1:
                return [ (face[0], face[1], closest0), (closest0, face[1], closest1) ]
            else:
                return [ (face[0], face[1], closest0) ]

    mesh = []
    for rgb,face in mesh0:
        for face2 in trimFace(face, polygon):
            mesh.append((rgb, face2))
    return mesh
    
if __name__ == '__main__':
    import cmath
    opts, args = getopt.getopt(sys.argv[1:], "e:UR:Uhdulw:P:o:Oc:LT:M:m:A:XHrf:na:D:t:s:S:x:y:z:Z:p:f:F:", 
                    ["help", "down", "up", "lower-left", "allow-repeats", "no-allow-repeats", "scale=", "config-file=",
                    "area=", 'align-x=', 'align-y=', 'optimization-time=', "pens=",
                    'input-dpi=', 'tolerance=', 'send=', 'send-speed=', 'work-z=', 'lift-delta-z=', 'safe-delta-z=',
                    'pen-down-speed=', 'pen-up-speed=', 'z-speed=', 'hpgl-out', 'no-hpgl-out', 'shading-threshold=',
                    'shading-angle=', 'shading-crosshatch', 'no-shading-crosshatch', 'shading-avoid-outline', 
                    'pause-at-start', 'no-pause-at-start', 'min-x=', 'max-x=', 'min-y=', 'max-y=',
                    'no-shading-avoid-outline', 'shading-darkest=', 'shading-lightest=', 'stroke-all', 'no-stroke-all', 'gcode-pause', 'dump-options', 'tab=', 'extract-color=', 'sort', 'no-sort', 'simulation', 'no-simulation', 'tool-offset=', 'overcut=', 
                    'boolean-shading-crosshatch=', 'boolean-sort=', 'tool-mode=', 'send-and-save=', 'direction=' ], )
    saveSCAD("circle.scad", inflatePolygon([ (30 * cmath.exp( 2 * 3.141519 * 1j * k / 40 ), 30 * cmath.exp( 2 * 3.141519 * 1j * (k+1) / 40 )) for k in range(40) ]))