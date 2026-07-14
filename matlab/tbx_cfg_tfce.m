function tfce_estimate = tbx_cfg_tfce
% SPM Configuration file for TFCE estimate
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________


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
if strcmp(spm('ver'),'SPM8')
  mask.filter  = {'image'};
else
  mask.filter  = {'image','mesh'};
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
% sequential stopping
% ---------------------------------------------------------------------
use_sequential_stopping        = cfg_menu;
use_sequential_stopping.tag    = 'use_sequential_stopping';
use_sequential_stopping.name   = 'Stop early if nothing can become significant';
use_sequential_stopping.labels = {'no','yes'};
use_sequential_stopping.values = {0 1};
use_sequential_stopping.val    = {0};
use_sequential_stopping.help   = {[...
'Stop the permutations as soon as it is certain that nothing in the image can become significant, instead of always running the full number of permutations.']
''
[...
'The largest value in the image is watched. No other element rests on fewer exceedances than it does, because an element with a smaller value is exceeded at least as often by the permutation maxima. Once the largest value has been exceeded often enough, and its corrected p-value is at least three standard errors above the largest alpha you asked about, no element in the image can still become significant and no corrected p-value can still move appreciably. The permutations are stopped there and the results are saved as usual.']
''
[...
'An image with nothing in it reaches that almost at once, because its largest value is an ordinary draw from the very distribution it is being compared against. An image with a real effect never reaches it, nor does an image whose corrected p-value sits anywhere near alpha, and both run the full number of permutations. So the permutations are only ever saved where they could not have changed the answer. At least 500 permutations are always run, which the tail approximations need.']
''
[...
'This is off by default: it changes how many permutations a given analysis runs, and analyses that were run with a fixed number of permutations are easier to compare with one another. Switch it on when you are screening many contrasts or designs and most of them are expected to be null.']
''
'Besag & Clifford (1991); Gandy (2009); Winkler et al. (2016), Faster permutation inference in brain imaging.'
};

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
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts n_perm};
conspec.help    = {''};

% ---------------------------------------------------------------------
% exchangeability blocks
% ---------------------------------------------------------------------
exch_blocks         = cfg_files;
exch_blocks.tag     = 'exch_blocks';
exch_blocks.name    = 'Exchangeability blocks';
exch_blocks.filter  = 'any';
exch_blocks.ufilter = '.*\.(txt|csv|mat)$';
exch_blocks.num     = [0 1];
exch_blocks.val     = {{}};
exch_blocks.help    = {[...
'Optionally define the exchangeability blocks explicitly. Provide a text file with one integer label per data point, in the order of the rows of the design matrix, or a mat-file that contains such a vector. ',...
'Data are then only permuted WITHIN blocks that share the same label. '],...
'',...
[...
'If left empty, the blocks are recognized automatically from the design: the subject effects of a paired t-test, a within-subject Anova or a flexible factorial define one block per subject, and all other designs are freely exchangeable. '],...
'',...
[...
'Use this option if your design is not recognized correctly, or to describe a block structure that is not part of the design matrix at all, such as scanning site or family membership. ',...
'Note that a one-sample t-test is permuted by sign-flipping and is therefore not affected by exchangeability blocks.']};

% ---------------------------------------------------------------------
% results Results Report
% ---------------------------------------------------------------------
tfce_estimate          = cfg_exbranch;
tfce_estimate.tag      = 'tfce_estimate';
tfce_estimate.name     = 'Estimate TFCE';
tfce_estimate.val      = {data nproc mask conspec nuisance_method use_sequential_stopping tbss exch_blocks};
tfce_estimate.help     = {''};
tfce_estimate.prog     = @tfce_estimate_stat;
