from vector import Vector
import math
import cmath
import triangulate

def reciprocalIntegral(p, a, b):
    # integrate 1/(z-p) from a to b
    return cmath.log((b-p)/(a-p))

def windingNumberSegments(p, segments, strict=False):
    segments = tuple((Vector(s[0]).toComplex(),Vector(s[1]).toComplex()) for s in segments)
    p = Vector(p).toComplex()
    for s in segments:
        if p in s:
            return not strict
    # we will integrate 1/(z-p) along the line segments; assume neither a nor b is 0
    total = sum(reciprocalIntegral(p, s[0], s[1]).imag for s in segments)
    return int(math.floor(0.5+total/(2*math.pi)))
    
def windingNumberLoop(p, loop, strict=False):
    loop = tuple(Vector(v).toComplex() for v in loop)
    p = Vector(p).toComplex()
    if p in loop:
        return not strict
    total = sum(reciprocalIntegral(p, loop[i-1], loop[i]).imag for i in range(len(loop)))
    return int(math.floor(0.5+total/(2*math.pi)))

if __name__ == '__main__':
    loops = 2
    points = 6
    segments = []

    for j in range(loops):
        r = (loops-j) * 40
        sign = -1 if j%2 else 1
        for i in range(points):
            angle1 = sign * 2 * math.pi * i / points
            angle2 = sign * 2 * math.pi * ((i+1)%points) / points
            segments.append( ( (r*math.cos(angle1),r*math.sin(angle1)), (r*math.cos(angle2),r*math.sin(angle2)) ) )

    print(windingNumberSegments(0, segments))
    print(windingNumberLoop(0, triangulate.lineSegmentsToLoop(segments)))
