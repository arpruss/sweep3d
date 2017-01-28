Python (2.7/3.x) code to sweep a changing cross-section along a smooth curve 
to generate either an STL mesh or an OpenSCAD source file. 

The OpenSCAD option works better if there is self-intersection---the sweep code isn't
smart enough to handle self-intersections.

The curve is divided into segments, and each segment generates a little straigh tube. You need
to make sure that the segments are long enough that the difference in curvature between the
ends doesn't make the two end-faces intersect. (To see how that could happen, increase the
size of the cross-section.)

There is also a simple vector and matrix library, and an ear-clipping 2D triangulator.

To try this out, just run sweep.py and you'll generate rings.stl and rings.scad, which
are Borromean rings with five-pointed-star cross-section.