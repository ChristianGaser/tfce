function tfce_estimate = tbx_cfg_tfce_estimate
% SPM Configuration file for TFCE estimate
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id$

rev = '$Rev$';

addpath(fileparts(which(mfilename)));

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file that contains the design specification.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

% ---------------------------------------------------------------------
% mask Select mask image to restrict analysis
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Select additional mask image';
mask.help    = {'Select an additional mask image to restrict analysis. As default the mask image in the analysis folder is used. Here you can select a mask image to additionally restrict the analysis to regions of interest (i.e. small volume correction).'};
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
titlestr.name    = 'Results Title';
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
n_perm.help    = {'With 1000 permutations the smallest possible p-value is 0.001 (n=1/p). A useful strategy is to start with 1000 permutations and continue to 5000-10000 only if p is small enough to be interesting and/or for the final analysis.'
                     ''
                     'If number of maximal possible permutations is smaller, then this number is used resulting in an exact permutation test.'
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
'Use 2D optimization (e.g. for TBSS data) with internal TFCE parameters H=2, E=1.']};

% ---------------------------------------------------------------------
% variance smoothing
% ---------------------------------------------------------------------
freedman_lane         = cfg_menu;
freedman_lane.tag     = 'freedman_lane';
freedman_lane.name    = 'Permutation method to deal with nuisance variables';
freedman_lane.labels = {'Draper-Stoneman','Freedman-Lane'};
freedman_lane.values  = {0 1};
freedman_lane.val     = {0};
freedman_lane.help    = {'A number of methods are available to obtain parameter estimates and construct a reference distribution in the presence of nuisance variables. Freedman-Lane permutation method can be optionally used if any nuisance variables exist. If no nuisance variables were found in the model then Draper-Stoneman method is automatically used.'
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
tfce_estimate.val      = {spmmat mask conspec freedman_lane tbss singlethreaded};
tfce_estimate.help     = {''};
tfce_estimate.prog     = @cg_tfce_estimate;
