function startup
% to add submodule: git submodule add https://github.com/flatironinstitute/FMM3D.git <destination_folder>
% might need to run: git submodule update --init --recursive

% compile rootfinder (c)
if 0 % add later
% if ~exist([pwd '/utils/bin/rootfinder_initial_guess_mex'],'file') % check 1st mex file...
  make_c
end
% compile kdtree (c++)
% if ~exist([pwd '/utils/bin/kdtree_ball_query.mex'],'file')
if 0 % add later
% if isempty(dir(fullfile(pwd,'/utils/bin/kdtree_ball_query.mex*')))
  make_cpp_kdtree
end
% compile close eval quadr (fortran)
if ~exist([pwd '/utils/bin/specialquad'],'file')
  disp('run makefile in /utils/f/ to build mex function for moments')
  if ismac
    if strcmp(computer,'MACI64') % intel
      !cp -f utils/f/make.inc.macosx_gcc utils/f/make.inc;    
    end
    if strcmp(computer,'MACA64')
      !cp -f utils/f/make.inc.macosx_arm64 utils/f/make.inc;    
    end
  elseif isunix
    !cp -f utils/f/make.inc.linux_gcc utils/f/make.inc;   
  end
  cd './utils/f';
  !make;
  cd '../..';
  !cp -f utils/f/specialquad.mex* utils/bin/
end
% compile FMM3D (fortran)
if isempty(dir(fullfile(pwd,'/utils/bin/fmm3d.mex*'))) 
  if ismac
    if strcmp(computer,'MACI64') % intel
      !cp -f utils/f/make.fmm3d.inc.macos.gnu utils/FMM3D/make.inc;    
    end
    if strcmp(computer,'MACA64')
      !cp -f utils/f/make.fmm3d.inc.macos_arm64.inc utils/FMM3D/make.inc;    
    end
  end
  cd './utils/FMM3D';
  !make clean;
  !make matlab;  
  cd '../..';
  !cp -f utils/FMM3D/matlab/fmm3d.mex* utils/bin/
end

