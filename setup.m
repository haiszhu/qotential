function setup
% Add the maintained qotential MATLAB paths.
mfilepath = fileparts(mfilename('fullpath'));
addpath(fullfile(mfilepath, 'utils'));
addpath(fullfile(mfilepath, 'utils', 'harmonic'));
addpath(fullfile(mfilepath, 'utils', 'bin'));
addpath(fullfile(mfilepath, 'matlab'));
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
end
