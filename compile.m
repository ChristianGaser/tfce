function compile(check)
% Compile the mex-files of the TFCE toolbox
% FORMAT compile
%        compile(1)   raise an error if a mex-file could not be built
%
% The optional check flag is used by the continuous integration, where a failed
% compilation has to abort the build instead of only printing a message.
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin < 1, check = 0; end

if strcmp(mexext,'mexmaci64') && verLessThan('matlab','9.2')
  warning('WARNING: Matlab version should be at least R2017a for compilation under Mac.');
end

if ispc
  % MSVC understands neither the gcc/clang flags nor pthreads. The mex-files use
  % Win32 threads there instead (see tfce_threads.h), and mex already optimizes.
  mexflag = ' -O ';
else
  mexflag = ' -O COPTIMFLAGS=''-O3 -fwrapv -DNDEBUG'' CFLAGS=''$CFLAGS -pthread -Wall -ansi -pedantic -Wextra'' ';
end

% name of the mex-file and a call that exercises it
tests = {
  'tfceMex_maxtree',       @() tfceMex_maxtree(rand(10,10,10),0.5,2,1)
  'tfceMex_maxtree_batch', @() tfceMex_maxtree_batch(rand(1000,4),0.5,2,1,[10 10 10])
  'tfceMex_resss',         @() tfceMex_resss(single(rand(100,10)),single(rand(100,2)),rand(10,2))
};

failed = {};

for i = 1:size(tests,1)
  name = tests{i,1};

  try
    eval(['mex ' mexflag ' ' name '.c'])
    tests{i,2}();
    fprintf('Compilation of %s successful\n', name);
  catch err
    fprintf('Compilation of %s not successful: %s\n', name, err.message);
    failed{end+1} = name; %#ok<AGROW>
  end
end

if check && ~isempty(failed)
  error('Could not build: %s', strjoin(failed, ', '));
end
