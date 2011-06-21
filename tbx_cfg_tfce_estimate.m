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
contrasts.name    = 'Contrast';
contrasts.help    = {
                     'Index of contrast.'
}';
contrasts.strtype = 'e';
contrasts.num     = [1 1];
% ---------------------------------------------------------------------
% number of permutations
% ---------------------------------------------------------------------
n_perm         = cfg_entry;
n_perm.tag     = 'n_perm';
n_perm.name    = 'Number of permutations';
n_perm.help    = {
                     'Number of permutations.'
                     ''
                     'If number of maximal possible permutations is smaller, then this number is used.'
}';
n_perm.strtype = 'e';
n_perm.num     = [1 Inf];
% ---------------------------------------------------------------------
% variance smoothing
% ---------------------------------------------------------------------
vFWHM         = cfg_entry;
vFWHM.tag     = 'vFWHM';
vFWHM.name    = 'Variance smoothing (for low DFs)';
vFWHM.help    = {
                     'Variance smoothing (for low DFs).'
                     ''
                     'For low DFs this option allows to smooth the variance.'
}';
vFWHM.strtype = 'e';
vFWHM.val     = {0};
vFWHM.num     = [0 Inf];
% ---------------------------------------------------------------------
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts n_perm vFWHM};
conspec.help    = {''};
% ---------------------------------------------------------------------
% results Results Report
% ---------------------------------------------------------------------
tfce_estimate          = cfg_exbranch;
tfce_estimate.tag      = 'tfce_estimate';
tfce_estimate.name     = 'Estimate TFCE';
tfce_estimate.val      = {spmmat conspec};
tfce_estimate.help     = {''};
tfce_estimate.prog     = @cg_run_tfce_estimate;
