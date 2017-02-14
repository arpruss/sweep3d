from surface import *
import svgpath.shader as shader
import svgpath.parser as parser
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
    bottomLeftV = Vector(bottomLeft)
    surface = inflateRaster(raster, thickness=thickness, roundness=roundness, iterations=iterations)
    mesh0 = surfaceToMesh(surface, center=False, twoSided=twoSided)

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
                #z = z0 + 0.5 * length * normDelta
                return (z.real, z.imag, 0.)
            else:
                print "notrim", start, stop, delta, rotate, state.bestLength
                return stop
    
        zeroCount = sum(1 for v in face if v.z == 0.)
        if zeroCount == 3 or zeroCount == 0:
            return [face]
        elif zeroCount == 2:
            if face[1].z != 0:
                face = (face[1], face[2], face[0])
            elif face[2].z != 0:
                face = (face[2], face[0], face[1])
            # now, the first vertex has non-zero z
            return [ (face[0], trimLine(face[0], face[1]), trimLine(face[0], face[2])) ]
        else: # zeroCount == 1
            if face[0].z == 0.:
                face = (face[1], face[2], face[0])
            elif face[1].z == 0.:
                face = (face[2], face[0], face[1])
            # now, the first two vertices have non-zero z
            closest0 = trimLine(face[0], face[2])
            closest1 = trimLine(face[1], face[2])
            if closest0 != closest1:
                return [ (face[0], face[1], closest0), (closest0, face[1], closest1) ]
            else:
                return [ (face[0], face[1], closest0) ]

    mesh = []
    for rgb,face in mesh0:
        face = tuple( Vector(v)*spacing+bottomLeftV for v in face )
        add = (face,)
        for v in face:
            if v.z == 0.:
                add = trimFace(face, polygon, shadeMode=shader.Shader.MODE_EVEN_ODD)
                break
        mesh += ((rgb,face) for face in add)
    return mesh
    
if __name__ == '__main__':
    import math
    import cmath
    shape = [ [30*cmath.exp(2*math.pi*1j*k/50),30*cmath.exp(2*math.pi*1j*(k+1)/50)] for k in range(50) ]
    saveSCAD("triangle.scad", inflatePolygon(shape,twoSided=True))