function ok = run_all
% Run the whole validation suite of the TFCE toolbox
% FORMAT ok = run_all
%
% Requires SPM on the path (for the surface checks) and the compiled mex-files.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fileparts(here));

tests = {
  'val_tfce_exactness'      % is the max-tree the exact TFCE integral?
  'val_gamma'               % is the Gamma approximation of the FWE null calibrated?
  'val_pareto'              % does the Pareto tail resolve uncorrected P below 1/n_perm?
  'val_sequential'          % does stopping early reach the same answer as the full run?
  'val_half_permutations'   % is the half-permutation shortcut lossless?
  'val_nuisance'            % are Draper-Stoneman / Freedman-Lane / Smith calibrated?
  'val_glm_fast'            % is the accelerated GLM identical to the unaccelerated one?
  'val_voxel_covariate'     % is the voxel-wise covariate path calibrated, and how slow?
};

failed = {};
t0 = tic;

for i = 1:numel(tests)
  try
    feval(tests{i});
  catch err
    fprintf('\n  *** %s raised an error: %s\n', tests{i}, err.message);
    failed{end+1} = tests{i}; %#ok<AGROW>
  end
end

fprintf('\n%s\n', repmat('=',1,74));
if isempty(failed)
  fprintf('  validation suite finished in %.0f s\n', toc(t0));
else
  fprintf('  the following scripts raised an error: %s\n', strjoin(failed, ', '));
end
fprintf('%s\n\n', repmat('=',1,74));

ok = isempty(failed);
