%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This matlab code was written as a part of my specialization project during the 
fall of 2015. The project is a part of my masters thesis and was written in 
cooperation with SINTEF ICT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are two gridding functions included:

compositePebiGrid.m
pebiGrid.m

Both functions creates a mrst grid structure ( http://www.sintef.no/contentassets/8af8db2e42614f7fb94fb0c68f5bc256/mrst-book-2015.pdf). The first creates a semi-structured grid, while the last creates a fully unstructured grid using distmesh. The last routine uses DistMesh: Per-Olof Persson and Gilbert
Strang, "A Simple Mesh Generator in MATLAB," SIAM Review Vol. 46 (2)
2004.

Example of how to run the code can be found in the script examples.m



