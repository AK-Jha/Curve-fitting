# Curve-fitting
This Matlab program uses the finite element method to fit a general curve.
The algorithm is designed to obtain an exact geometrical representation of the boundary curves (or the geometry itself) for use in finite element analysis to perform engineering calculations in complex domain boundaries.
After further optimization, it can also be used for general purposes such as graphic design or computer graphics.

One of the best outcomes of this technique is the automatic care of the curvature of the curve: where the curvature is high, the segment size is small, and where the curvature is flat, the segment size is large.
This determines the element size in the mesh of the two- or three-dimensional geometries represented by these curves, which is small at complex curve boundaries and large in straight regions.
When these fitted points (along with the approximation of fit) are used to create the mesh, the error in geometry stays within the desired error norm.
This gives high accuracy, lowers computation costs, and gets rid of the usual practice of remeshing geometry near curved edges to get a better representation.
