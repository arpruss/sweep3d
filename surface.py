from vector import *
from exportmesh import *
from PIL import Image
import itertools
import os.path

def surfaceToMesh(data, center=False, twoSided=False, zClip=None, xScale=1., yScale=1., zScale=1., tolerance=0., color=None):
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
        offsetX = 0
        offsetY = 0
        
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

def inflateRaster(raster, thickness=10., roundness=1., iterations=None):
    """
    raster is a boolean matrix.
    
    roundness varies from 0.0 for maximally steep sides to 1.0 for a very gradual profile.
    
    Here's a way to visualize how inflateImage() works. The white or transparent areas 
    of the image are cold, clamped at temperature 0. Above the image, there is a layer of
    insulation, and above that there is a flat heater at a fixed positive temperature. So, 
    the clamped areas stay at the fixed temperature, but away from the clamped areas, the image
    heats up. How the heat profile looks depends on how effective the layer of insulation is.
    The effectiveness of the layer of insulation is measured by the roundness parameter.
    When this parameter is close to zero, the unclamped areas all get to close to the heater
    temperature, resulting in a very sharp heat gradient near the image boundaries, and flatness
    further away from the boundaries. When the parameter is close to 1, the insulation is more
    effective, and the temperature rise away from the boundaries is more gradual. The result is
    a smoother change in the gradient.
    
    Finally, the temperatures are scaled to be between 0 and thickness to generate the inflated
    height map. 
    
    If roundness > 1., things blow up due to uncontrolled heating.
    """
    width = len(raster)
    height = len(raster[0])
        
    if iterations == None:
        iterations = 60 * max(width,height)
        
    data = [ [0. for y in range(height)] for x in range(width) ]
    
    for iter in range(iterations):
        newData = [ [0. for y in range(height)] for x in range(width) ]
        for x in range(width):
            for y in range(height):
                def z(dx,dy):
                    x1 = x+dx
                    y1 = y+dy
                    if x1 < 0 or x1 >= width or y1 < 0 or y1 >= height:
                        return 0.
                    else:
                        return data[x1][y1]
                        
                if raster[x][y]:
                    newData[x][y] = 1.0+roundness*(z(-1,0)+z(1,0)+z(0,-1)+z(0,1)+0.7*(z(-1,-1)+z(1,1)+z(-1,1)+z(1,-1)))/(4+0.7*4)

        data = newData
        
    maxZ = max(max(col) for col in data)
    
    return [ [datum / maxZ * thickness for datum in col] for col in data ]


def inflateImage(image, thickness=10., roundness=1., iterations=None):
    def inside(x,y):
        rgb = image.getpixel((x,image.size[1]-1-y))
        if len(rgb) > 3 and rgb[3] == 0:
            return False
        return rgb[0:3] != (255,255,255)

    raster = [ [ inside(x,y) for y in range(image.size[1]) ] for x in range(image.size[0]) ]
    
    return inflateRaster(image, thickness=thickness, roundness=roundness, iterations=itierations)
        
if __name__ == '__main__':
    inPath = sys.argv[1]
    outPath = os.path.splitext(inPath)[0] + ".scad"
    baseName = os.path.splitext(os.path.basename(outPath))[0]
    
    thickness = 10.
    roundness = 1.
    if len(sys.argv)>2:
        thickness = float(sys.argv[2])
    if len(sys.argv)>3:
        roundness = float(sys.argv[3])

    image = Image.open(inPath).convert('RGBA')
    
    print("Inflating...")
    data = inflateImage(image,thickness=thickness,roundness=roundness)
    
    scadModule = toSCADModule(surfaceToMesh(data, twoSided=False, center=True), baseName+"_raw")
    scadModule += """

module %s() {
     render(convexity=2)
     translate([0,0,-%f])
     intersection() {
        %s_raw();
        translate([-%d/2,-%d/2,%f]) cube([%d,%d,%f]);
     }
}

%s();
""" % (baseName, thickness / 20., baseName, image.size[0], image.size[1], thickness / 20., image.size[0]+2, image.size[1]+2, thickness+1., baseName)
    
    print("Saving "+outPath)
    with open(outPath, "w") as f: f.write(scadModule)
