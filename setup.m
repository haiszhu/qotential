function setup
% all path...
disp('might need to run: git submodule update --init --recursive')
disp('then startup.m')
mfilepath=fileparts(mfilename('fullpath'));
addpath([mfilepath, '/kernels']);
addpath([mfilepath, '/utils']);
addpath([mfilepath, '/utils/harmonics']);
addpath([mfilepath, '/utils/bin']);

