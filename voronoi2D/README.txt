%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This matlab code was written as a part of my specialization project during the 
fall of 2015. The project is a part of my masters thesis and was written in 
cooperation with SINTEF ICT.
Runar Lie Berge                                      (runarlb@stud.ntnu.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are two functions included for generating unstructured in MRST 
(http://www.sintef.no/contentassets/8af8db2e42614f7fb94fb0c68f5bc256/mrst-book-2015.pdf).

Functions:
    compositePebiGrid.m
    pebiGrid.m
    example.m

compositePebiGrid.m creates a semi-structured grid, by inserting voronoi
seeds around wells and fractures.

pebiGrid.m creates a fully unstructured grid. It uses the software DistMesh.
DistMesh is a software for creating unstructured delaunay triangulations:
Per-Olof Persson and Gilbert Strang, "A Simple Mesh Generator in MATLAB,"
SIAM Review Vol. 46 (2) 2004.

Example of how to run the code can be found in the script examples.m.
