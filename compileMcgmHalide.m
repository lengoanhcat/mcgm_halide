%% (re)Compile mcgmOpticalFlow.mexa64 with Halide

% Add the path to mex_halide.m.
% setenv('HALIDE_PATH','/usr/local/halide'); %% previously used
setenv('HALIDE_PATH','/home/lpzal1/halide_autoscheduler.build/');
addpath(fullfile(getenv('HALIDE_PATH'), 'tools'));

% $$$ setenv('CPP_INCLUDE_PATH',['/usr/local/matlab/R2016a/extern/' ...
% $$$                     'include']);
% $$$ setenv('LIBRARY_PATH','/usr/local/matlab/R2016a/bin/glnxa64/');

% Build the mex library from the blur generator.
setenv('HL_DEBUG_CODEGEN','2');
setenv('HL_NUM_THREADS','12');
% $$$ setenv('HL_TRACE','1');
% $$$ setenv('HL_TRACE_FILE','./mcgmOpticalFlow.trace');
% mex_halide('./mcgmOpticalFlow.cpp','-e html');
%% compile v01 version of motion model
% mex_halide('./mcgmOpticalFlow_v01.cpp')
%% compile v01 version of stereo model
% mex_halide('./mcgmOpticalFlow_v01_stereo.cpp');
%% compile v02 version of motion model
% mex_halide('./mcgmOpticalFlow_v02.cpp')
%% compile v02 version of motion model with autoschedule
mex_halide('./mcgmOpticalFlow_v02_autoscheduled.cpp')
%% compile v02 version of stereo model
% mex_halide('./mcgmOpticalFlow_v02_stereo.cpp'); % compile v0.3 version of motion model
% mex_halide('./mcgmOpticalFlow.cpp','-e bitcode');

%% GPU
% mex_halide('./mcgmOpticalFlow_gpu.cpp');
