from vector import *
from exportmesh import *
from PIL import Image
import itertools

def surfaceToMesh(data, center=False, twoSided=False, zClip=None, xScale=1., yScale=1., zScale=1., tolerance=1e-6, color=None):
    # clipping is done before z-scaling
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
            elif zClip is None:
                return data[x][y]
            else:
                return min(data[x][y],zClip)

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
    
    # this code could be made way more efficient
    
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

def inflateImage(image, thickness=10., roundness=1., iterations=None):
    """
    roundness varies from 0.0 for maximally steep sides to 1.0 for a very smooth profile 
    """
    def inside(x,y):
        if x < 0 or x >= image.size[0] or y < 0 or y >= image.size[1]:
            return false
        rgb = image.getpixel((x,image.size[1]-1-y))
        if len(rgb) > 3 and rgb[3] == 0:
            return False
        return rgb[0:3] != (255,255,255)
        
    width = image.size[0]
    height = image.size[1]
    
    if iterations == None:
        iterations = 60 * max(width,height)
        
    data = [ [0. for y in range(height)] for x in range(width) ]
    
    for iter in range(iterations):
        newData = [ [0. for y in range(height)] for x in range(width) ]
        for x in range(width):
            for y in range(height):
                def z(dx,dy):
                    if x+dx < 0 or x+dx >= width or y+dy < 0 or y+dy >= height:
                        return 0.
                    else:
                        return data[x+dx][y+dy]
                        
                if inside(x,y):
                    newData[x][y] = 1.0+roundness*(z(-1,0)+z(1,0)+z(0,-1)+z(0,1)+0.7*(z(-1,-1)+z(1,1)+z(-1,1)+z(1,-1)))/(4+0.7*4)
        data = newData
        
    maxZ = max(max(col) for col in data)
    
    return [ [datum / maxZ * thickness for datum in col] for col in data ]
        
if __name__ == '__main__':
    data = [ [0.1,0.2], [0.3,0.4] ]
    saveSTL("quad.stl", surfaceToMesh(data, twoSided=True))
    
    image = Image.open('smallheart.png').convert('RGBA')
    
    print("Inflating...")
    data = inflateImage(image,thickness=10,roundness=0.98,iterations=3000)
    
    scadModule = toSCADModule(surfaceToMesh(data, twoSided=False), "smallHeart")
    scadModule += """

    render(convexity=2)
    translate([0,0,-0.5])
    intersection() {
        smallHeart();
        translate([0,0,.5]) cube([200,200,200]);
    }
"""
    
    print("Saving smallheart.scad")
    with open("smallheart.scad", "w") as f: f.write(scadModule)
    