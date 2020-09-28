# VTK_python
Useful VTK functions implemented in Python </br>
This script contains the following functions:
- read vtu() : read the unstructured data
- read vtp() : read Polygonal Data
- get_normals() : return the points normal of the VTK file
- get_curvature_mean(): return the mean curvature of the points.
- get_curvature_gaussian() :return the gaussian curvature of the points.
- get_neighbors_around_minima() : take a point as first argument and returns the nearest K neighbors to this point (K can be determined inside the function)
- visualize_data() : visualize the vtk model (VTK has a nice software for that called "Paraview")
- get_minimum() : if vtk data has a scalar value (for example safety factor in engineering application), then it is needed to find the points with the lowest SF. this function returns the points which has the lowest value. 
- write_vtu_ascii() : write the data into vtu file in an Ascii mode.
- write_vtp() : write the data into vtp file.
