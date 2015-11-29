BEFORE YOU START:
to run some of the functions you need the tweeked version of the 3rd party module distmesh. This can be found in the distmesh/ folder. To use it you have to let mrst know where to find it. Run the following lines in matlab (Remark: you have to run it from the same directory you found this readme file)

path = fullfile(pwd,'distmesh/')
mrstPath('reregister','distmesh',path)

You can copy these lines into your startup.m file, so they run whenever you start matlab. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


There are two gridding functions included;

compositeGridPEBI.m
compositeGridPEBIdistmesh.m

The first creates a semi-structured grid, while the last creates a fully unstructured grid using distmesh.


Examples using the code can be found in the folder testing (NOTE, some of the testing scripts may be outdated).

