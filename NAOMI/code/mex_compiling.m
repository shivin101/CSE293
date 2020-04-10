%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEX SCRIPT
% 
% This script will compike all the mex files needed for the NAOMI simulator
%
% 2017 - Adam Charles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ./MEX/
fprintf('Compiling array_SubMod.cpp and array_SubModTest.cpp for faster array substitution...\n')
mex -largeArrayDims array_SubMod.cpp
mex -largeArrayDims array_SubModTest.cpp
fprintf('Compiling array_SubSub.cpp and array_SubSubTest.cpp for faster array summing...\n')
mex -largeArrayDims array_SubSub.cpp
mex -largeArrayDims array_SubSubTest.cpp
fprintf('Compiling dendrite_dijkstra_cpp.cpp for faster dendrite growth...\n')
mex -largeArrayDims dendrite_dijkstra_cpp.cpp
fprintf('Compiling dendrite_randomwalk_cpp.cpp for faster dendrite growth...\n')
mex -largeArrayDims dendrite_randomwalk_cpp.cpp
cd ..
fprintf('done.\n')