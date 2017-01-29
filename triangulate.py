from vector import *
from random import shuffle

# triangulation algorithm based on https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
# minus all the double-linked list stuff that would be great but I am not bothering with it

def cross_z(a,b):
    # z-component of cross product
    return a.x * b.y - a.y * b.x

def pointInside(p,a,b,c):
    # checks if point is strictly inside the triangle a,b,c
    c1 = cross_z(p-a, b-a)
    c2 = cross_z(p-b, c-b)
    if c1 * c2 <= 0:
        return False
    c3 = cross_z(p-c, a-c)
    return c1 * c3 > 0 and c2 * c3 > 0
    
def triangulate(polygon):
    # assume input polygon is counterclockwise
    n = len(polygon)

    if n < 3:
        raise Exception
    
    # two efficient special cases
    if n == 3:
        return [(0,1,2)]
    elif n == 4:
        return [(0,1,2),(2,3,0)]
    
    triangles = []
    polygon = [Vector(v) for v in polygon]
    
    index = list(range(n))
    
    def isReflex(i):
        return cross_z(polygon[i]-polygon[i-1], polygon[(i+1) % n]-polygon[i]) < 0
        
    reflex = [isReflex(i) for i in range(n)]
    
    def isEar(i):
        if reflex[i]:
            return False
        a,b,c = polygon[i-1],polygon[i],polygon[(i+1) % n]
        j = i+2
        while j % n != (i-1) % n:
            if reflex[j % n] and pointInside(polygon[j % n],a,b,c):
                return False
            j += 1
        return True
        
    ear = [isEar(i) for i in range(n)]
    
    while n >= 3:
        foundEar = False
        # the backwards search makes the deletions faster
        # if we had a double-indexed list, we wouldn't have this to worry about
        for i in range(n-1,-1,-1):
            if ear[i]:
                triangles.append((index[i-1],index[i],index[(i+1) % n]))
                # TODO: less deleting!
                del polygon[i]
                del reflex[i]
                del ear[i]
                del index[i]
                n -= 1
                # it's tempting to optimize here for the case where n==3, but that would
                # probably be counterproductive, as the n==3 test would run O(n) times
                if reflex[i-1]:
                    reflex[i-1] = isReflex(i-1)
                if reflex[i % n]:
                    reflex[i % n] = isReflex(i % n)
                ear[i-1] = isEar(i-1)
                ear[i % n] = isEar(i % n)
                foundEar = True
                break
        assert foundEar
    return triangles

def polygonsToSVG(vertices, polys):
    vertices = tuple(Vector(v) for v in vertices)
    minX = min(v.x for v in vertices)
    minY = min(v.y for v in vertices)
    maxX = max(v.x for v in vertices)
    maxY = max(v.y for v in vertices)

    svgArray = []
    svgArray.append('<?xml version="1.0" standalone="no"?>')
    svgArray.append('<svg width="%fmm" height="%fmm" viewBox="0 0 %f %f" xmlns="http://www.w3.org/2000/svg" version="1.1">'%(maxX-minX,maxY-minY,maxX-minX,maxY-minY))
    for p in polys:
        path = '<path stroke="black" stroke-width="0.25px" fill="yellow" d="'
        for i in range(len(p)+1):
            path += 'L' if i else 'M'
            path += '%.6f %.6f ' % ( vertices[p[i % len(p)]] - (minX,minY) )
        path += '/>'
        svgArray.append(path)
    svgArray.append('</svg>')
    return '\n'.join(svgArray)
    
if __name__ == '__main__':
    import cmath
    import math
    import sys
    import time
    if len(sys.argv) >= 2:
        m = int(sys.argv[1])
    else:
        m = 16
    polygon = [ cmath.exp(2j * math.pi * k / m) * (10 if k%2 else 18) for k in range(m) ]
    t = time.time()
    tr = triangulate(polygon)
    t = time.time() - t
    sys.stderr.write("Time %.4fs\n" % t)
    #sys.stderr.write(str(tr)+"\n")
    print polygonsToSVG(polygon, tr)
    