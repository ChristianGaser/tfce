tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom, n_threads)
% Apply exact TFCE to many maps at once, one map per thread
% FORMAT tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom, n_threads)
% T         - N x B matrix, one T map per column
% E         - TFCE parameter for extent
% H         - TFCE parameter for height
% calc_neg  - also calc neg. TFCE values (default)
% geom      - [nx ny nz] for volume data, or an F x 3 face list for surfaces
% n_threads - number of worker threads (default: one per map)
%
% tfce      - N x B matrix of TFCE maps, one per column
%
% Intended for permutation testing. Parallelising across permutations instead of
% within a single map keeps every thread on independent data, so there is no
% locking on the hot path. Results are identical to calling tfceMex_maxtree once
% per map.
%
% Note that T and the result each need N*B*8 bytes, so B should be chosen in the
% order of the number of cores rather than as large as possible.
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

disp('Compiling tfceMex_maxtree_batch.c')

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

eval([mexcmd ' ' mexflag ' tfceMex_maxtree_batch.c'])

cd(p_path);

if nargin > 5
  tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom, n_threads);
else
  tfce = tfceMex_maxtree_batch(T, E, H, calc_neg, geom);
end

return
