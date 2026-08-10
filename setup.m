function setup
% all path...
disp('might need to run: git submodule update --init --recursive')
disp('then startup.m')
mfilepath=fileparts(mfilename('fullpath'));
addpath([mfilepath, '/kernels']);
addpath([mfilepath, '/utils']);
addpath([mfilepath, '/utils/harmonics']);
addpath([mfilepath, '/utils/bin']);
addpath([mfilepath, '/matlab']);
addpath(fullfile(mfilepath, 'external', 'LineQuaaadrature', 'matlab'));
addpath(fullfile(mfilepath, 'external', 'LineQuaaadrature', 'utils'));

fmm3d = getenv('FMM3D_DIR');
if isempty(fmm3d)
  fmm3d = fullfile(getenv('HOME'), 'git', 'FMM3D');
end
if isfolder(fullfile(fmm3d, 'matlab'))
  addpath(fullfile(fmm3d, 'matlab'));
else
  warning('setup:noFMM3D', ...
      'FMM3D not found at %s -- Lap3dSLPfmm and Lap3dDLPfmm will fail.', ...
      fmm3d);
end
