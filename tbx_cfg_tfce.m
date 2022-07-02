function tfce_estimate = tbx_cfg_tfce
% SPM Configuration file for TFCE estimate
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: tbx_cfg_tfce.m 215 2021-04-15 17:48:39Z gaser $

rev = '$Rev: 215 $';

addpath(fileparts(which(mfilename)));

% try to estimate number of processor cores
try
  if strcmpi(spm_check_version,'octave')
    numcores = nproc;
  else
    numcores = feature('numcores');
  end

  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
% or for Octave
if numcores == 1 || isdeployed || strcmpi(spm_check_version,'octave'), numcores = 0; end

%_______________________________________________________________________
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.hidden  = numcores <= 1 || isdeployed;
nproc.help    = {
    'In order to use multi-threading the TFCE job with multiple SPM.mat files can be split into separate processes that run in the background. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process might need a large amount of RAM, which should be considered to choose the appropriate number of processes.'
    ''
    'Please further note that additional modules in the batch can now be used because the processes are checked every minute.'
  };

% ---------------------------------------------------------------------
% data Select SPM.mat
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Select (v)SPM.mat';
data.help    = {'Select the (v)SPM.mat files that contain the design specification from a previous (parametric) estimation, where all required contrasts are already specified.'};
data.filter  = 'mat';
data.ufilter = 'SPM\.mat$';
data.num     = [1 Inf];

% ---------------------------------------------------------------------
% mask Select mask to restrict analysis
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Select additional mask';
mask.help    = {'Select an additional mask image or surface to restrict your analysis. As default the mask in the analysis folder is used. Here you can select a mask to additionally restrict the analysis to regions of interest (i.e. small volume/surface correction).'};
if strcmp(spm('ver'),'SPM12')
  mask.filter  = {'image','mesh'};
else
  mask.filter  = {'image'};
end
mask.val     = {''};
mask.ufilter = '.*';
mask.num     = [0 1];

% ---------------------------------------------------------------------
% titlestr Results Title
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'Results title';
titlestr.help    = {'Heading on results page - determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];

% ---------------------------------------------------------------------
% contrasts Contrast
% ---------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast index';
contrasts.help    = {'Index(es) of contrast according to the contrast manager.'
                     ''
                     'Each contrast in SPM is indicated by a sequential number that is displayed in the first column of the contrast manager.'
                     ''
                     'You can enter one or more contrasts. If only one number is entered, and this number is "Inf", you can select one or more contrasts interactively using the contrast manager.'
                     ''
                     'Do not define here the contrast itself. This should be done in the contrast manager, that is automatically called if "Inf" is kept as entry.'
}';
contrasts.strtype = 'e';
contrasts.val     = {Inf};
contrasts.num     = [1 Inf];

% ---------------------------------------------------------------------
% number of permutations
% ---------------------------------------------------------------------
n_perm         = cfg_entry;
n_perm.tag     = 'n_perm';
n_perm.name    = 'Number of permutations';
n_perm.help    = {'In order to obtain reliable estimates you need about 5000-10000 permutations.'
                  ''
                  'There is also an option to interrrupt the permutation process and to save the results at this step to take a first look on your results.'
                  ''
                  'If the number of maximal possible permutations is smaller, then this number is used resulting in an exact permutation test.'
                  ''
                  'Please note, that a tail approximation is finally used to estimate the corrected p-values. Thus, there is no dependency anymore between the lowest achievable p-value and the number of permutations as in previous versions.'
}';
n_perm.strtype = 'e';
n_perm.val     = {5000};
n_perm.num     = [1 Inf];

% ---------------------------------------------------------------------
% two-dimensional processing
% ---------------------------------------------------------------------
tbss    = cfg_menu;
tbss.tag = 'tbss';
tbss.name = 'TBSS data';
tbss.labels = {'yes','no'};
tbss.values = {1 0};
tbss.val  = {0};
tbss.help = {[...
'Use 2D optimization (e.g. for TBSS DTI data) with internal TFCE parameters H=2, E=1.']};

% ---------------------------------------------------------------------
% method to deal with nuisance variables
% ---------------------------------------------------------------------
nuisance_method         = cfg_menu;
nuisance_method.tag     = 'nuisance_method';
nuisance_method.name    = 'Permutation method to deal with nuisance variables';
nuisance_method.labels = {'Draper-Stoneman','Smith','Freedman-Lane (experimental)'};
nuisance_method.values  = {0 2 1};
nuisance_method.val     = {2};
nuisance_method.help    = {'A number of methods are available to obtain parameter estimates and construct a reference distribution in the presence of nuisance variables. Smith permutation method is used if any nuisance variables exist and is selected by default. If no nuisance variables were found in the model then Draper-Stoneman method is automatically used. '
''
'Freedman-Lane is another permutation method to deal with nuisance parameters. However, behaviour of that method was found to be strange under some circumstances and you have to apply this method very carefully. '
''
'It is only necessary to change the permutation method if a large discrepancy between parametric and non-parametric statistic was found, which is indicated at the Matlab command line. '
}';

% ---------------------------------------------------------------------
% weighting of cluster-size
% ---------------------------------------------------------------------
E_weight         = cfg_menu;
E_weight.tag     = 'E_weight';
E_weight.name    = 'Weighting of cluster size';
E_weight.labels = {'More weighting of focal effects (E=0.5)','More weighting of larger clusters (E=0.6)','Very large weighting of larger clusters (E=0.7)'};
E_weight.values  = {0.5 0.6 0.7};
E_weight.val     = {0.5};
E_weight.help    = {'The idea of the TFCE approach is to combine focal effects with large voxel height as well as broad effects. The weighting of these effects is defined using the parameters E (extent) and H (height). Smith and Nichols (Neuroimage 2009) empirically estimated E=0.5 and H=2 for volume data to provide good statistical power. However, the empirically derived values found to be very sensitive for local effects, but not for broader effects. Thus, you can try to change the weighting if you expect more broader effects or even very broad effects by changing the weighting parameter E.'
                     ''
'Please note that for surfaces and TBSS data the weightings are different from that of 3D volume data with E=1 and H=2 and will be fixed and not changed by this setting.'
}';

% ---------------------------------------------------------------------
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts n_perm};
conspec.help    = {''};

% ---------------------------------------------------------------------
% multithreading
% ---------------------------------------------------------------------
singlethreaded    = cfg_menu;
singlethreaded.tag = 'singlethreaded';
singlethreaded.name = 'Use multi-threading to speed up calculations';
singlethreaded.labels = {'yes','no'};
singlethreaded.values = {0 1};
if ispc
  singlethreaded.val  = {1};
else
  singlethreaded.val  = {0};
end
singlethreaded.help = {[...
'Multithreading can be used to distribute calculations to multiple processors. ',...
'This will minimize calculation time by a large amount, but makes trouble on Windows machines, where it is deselected by default. ']};

% ---------------------------------------------------------------------
% results Results Report
% ---------------------------------------------------------------------
tfce_estimate          = cfg_exbranch;
tfce_estimate.tag      = 'tfce_estimate';
tfce_estimate.name     = 'Estimate TFCE';
tfce_estimate.val      = {data nproc mask conspec nuisance_method tbss E_weight singlethreaded};
tfce_estimate.help     = {''};
tfce_estimate.prog     = @tfce_estimate_stat;
