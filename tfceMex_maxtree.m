tfce = tfceMex_maxtree(t, E, H, calc_neg, faces)
% Apply threshold-free cluster enhancement (TFCE), computed exactly
% FORMAT tfce = tfceMex_maxtree(t, E, H, calc_neg, faces)
% t         - T map (3D volume, or Nx1 vector for surface data)
% E         - TFCE parameter for extent
% H         - TFCE parameter for height
% calc_neg  - also calc neg. TFCE values (default)
% faces     - surface faces; omit or leave empty for volume data
%
% The last argument selects the neighbourhood: if faces are omitted or empty,
% t is treated as a 3D volume with 26-connectivity, otherwise as surface data
% on the mesh defined by these faces.
%
% The TFCE integral is evaluated exactly using a max-tree (component tree) built
% with union-find, so there is no step size dh and no discretisation error. Run
% time is independent of any precision parameter.
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev$';

disp('Compiling tfceMex_maxtree.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);

if strcmpi(spm_check_version,'octave')
  mexcmd = 'mkoctfile --mex';
  mexflag=' -O -DOCTAVE ';
else
  mexcmd = 'mex';
  mexflag=' -O -largeArrayDims COPTIMFLAGS=''-O3 -fwrapv -DNDEBUG'' CFLAGS=''$CFLAGS -pthread -Wall -ansi -pedantic -Wextra'' ';
end

eval([mexcmd ' ' mexflag ' tfceMex_maxtree.c'])

cd(p_path);

if nargin > 4
  tfce = tfceMex_maxtree(t, E, H, calc_neg, faces);
else
  tfce = tfceMex_maxtree(t, E, H, calc_neg);
end

return
