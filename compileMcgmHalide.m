%% (re)Compile mcgmOpticalFlow.mexa64 with Halide

% Add the path to mex_halide.m.
setenv('HALIDE_PATH','/usr/local/halide');
addpath(fullfile(getenv('HALIDE_PATH'), 'tools'));

% $$$ setenv('CPP_INCLUDE_PATH',['/usr/local/matlab/R2016a/extern/' ...
% $$$                     'include']);
% $$$ setenv('LIBRARY_PATH','/usr/local/matlab/R2016a/bin/glnxa64/');

% Build the mex library from the blur generator.
setenv('HL_DEBUG_CODEGEN','1');
setenv('HL_NUM_THREADS','12');
% $$$ setenv('HL_TRACE','1');
% $$$ setenv('HL_TRACE_FILE','./mcgmOpticalFlow.trace');
% mex_halide('./mcgmOpticalFlow.cpp','-e html');
mex_halide('./mcgmOpticalFlow.cpp');
% mex_halide('./mcgmOpticalFlow.cpp','-e bitcode');
