function tfce_estimate_stat(job)
% main TFCE function for estimating TFCE statistics
%
% FORMAT tfce_estimate_stat(job)
% job        - job from tbx_cfg_tfce
%
% job fields:
%   data
%
%   nproc
%     0 - no multi-threading
%     x - number of processors for multi-threading
%
%   nuisance_method 
%     method to deal with nuisance variables
%     0 - Draper-Stoneman
%     1 - Freedman-Lane
%     2 - Smith
% 
%   single-threaded
%     0 - multi-threaded TFCE estimation
%     1 - single-threaded TFCE estimation
%
%   conspec.contrasts
%     Inf - interactive selection of contrast
%     x   - index of contrast(s)
%
%   mask
%     mask for restricting TFCE estimation (SVC)
%
%   tbss
%     TBSS weighting
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

global old_method_stat

% disable parallel processing for only one SPM.mat file
if numel(job.data) == 1
  job.nproc = 0;
end

% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && exist('cat_parallelize','file')
  cat_parallelize(job,mfilename,'data');
  return
% or run though all data step-by-step
elseif numel(job.data)>1
  data = job.data;
  for i=1:numel(job.data)
    job = rmfield(job,'data');
    job.process_index = i;
    job.data{1} = data{i};
    tfce_estimate_stat(job)
  end
  return
end

% Use tail approximation from the Gamma distribution for corrected P-values 
use_gamma_tail_approximation = true;
tail_approximation_wo_unpermuted_data = false;

% single-threaded?
singlethreaded = job.singlethreaded;

% convert to z-statistic
convert_to_z = false;

% option to stop estimation if no FWE-corrected result was found after defined number of permutations
stop_if_no_FWEeffects_found = Inf;

% use (too liberal) method for estimating maximum statistic from old release 
% r184 for compatibility purposes only that was estimating max/min statistics
% only inside pos./neg. effects and not both
if isfield(job,'old_method_stat')
  old_method_stat = job.old_method_stat;
elseif exist('old_method_stat','var') && isempty(old_method_stat)
  old_method_stat = false;
end

if old_method_stat
  fprintf('Use old method from r184 to estimate maximum statistics which might be too liberal.\n');
end

% save null distribution
save_null_distribution = true;

% method to deal with nuisance variables
% 0 - Draper-Stoneman
% 1 - Freedman-Lane
% 2 - Smith
nuisance_method = job.nuisance_method;

% display permuted design matrix (otherwise show t distribution)
show_permuted_designmatrix = true;

% allow to test permutations without analyzing data
% test mode is automatically chosen if only a SPM.mat is available
% without any data
test_mode = false;

% define stepsize for tfce
n_steps_tfce = 100;
    
% colors and alpha levels
col   = [0.25 0 0; 1 0 0; 1 0.75 0];
alpha = [0.05      0.01   0.001];
    
% give same results each time
if exist('rng','file') == 2
  rng('default')
  rng(0)
else
  rand('state',0);
end

% tolerance for comparing real numbers
tol = 1e-4;
    
% check spm version
if ~strcmp(spm('ver'),'SPM12')
  error('Please use any TFCE version < r215 that still supports SPM8.')
end

% indicate if voxel-wise covariates were modeled
[pt,nm] = fileparts(job.data{1});
is_vSPM = strcmp(nm,'vSPM');

load(job.data{1});
cwd = fileparts(job.data{1});

%-Check that model has been estimated
if ~isfield(SPM, 'xVol')
  str = { 'This model has not been estimated.';...
              'Would you like to estimate it now?'};
  if spm_input(str,1,'bd','yes|no',[1,0],1)
    cd(cwd)
    SPM = spm_spm(SPM);
  else
    return
  end
end
    
Ic0 = job.conspec.contrasts;

% if just one contrast is defined use this and skip interactive selection
if (numel(Ic0) == 1) && ~isfinite(Ic0) && isfield(SPM,'xCon') && (numel(SPM.xCon) == 1)
  Ic0 = 1;
end

% check whether contrast are defined
if (numel(Ic0) == 1)
  if ~isfinite(Ic0) || ~isfield(SPM,'xCon') || (isfield(SPM,'xCon') && isempty(SPM.xCon))
    [Ic0,xCon] = spm_conman(SPM,'T&F',Inf,...
          '  Select contrast(s)...',' ',1);
    SPM.xCon = xCon;
  end
end

% for default SPM.mat results has to be called first
if isempty(SPM.xCon(Ic0(1)).eidf) && ~is_vSPM
  fprintf('You have to call results first.\n');
  cat_spm_results_ui('Setup',SPM);
  load(job.data{1});
end

% check that no temporal filter was used
if isstruct(SPM.xX.K)
  fprintf('ERROR: No first level analysis with temporal correlations allowed.\n');
  return
end

% get some parameters from SPM
xX     = SPM.xX;
VY     = SPM.xY.VY;
n_data = size(xX.X,1);

% sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design
xX = correct_xX(xX);

% check whether voxel-wise covariate is defined
voxel_covariate = 0;
for i = 1:numel(SPM.xC)
  if isfield(SPM.xC(i),'P') && ~isempty(SPM.xC(i).P)
    voxel_covariate = i;
    fprintf('Voxel-wise covariate found.\n');
  end
end

% find exchangeability block labels for longitudinal designs (paired t-test, flexible factorial)
repeated_anova = ~isempty(xX.iB);
if repeated_anova
  [rw,cl] = find(xX.I == length(xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_subj)
  exch_block_labels = xX.I(:,cl(1));     % column from above contains the subject factor

  % check that labels are defined for each block
  for i=1:n_data
    groupListed(exch_block_labels(i)) = true;
  end
  
  for i = 1:max(exch_block_labels)
    if ~groupListed(i) 
      fprintf('\nError: block %d must be assigned to at least one design row in the blocks file.\n\n',i);
      return
    end
  end
  
  % check whether ate least two time points are available
  n_tp = sum(xX.X(:,xX.iB));
  if any(n_tp == 1)
    fprintf('\nError: You need at least two time points for a longitudinal design for every subjects.\n\n');
    return
  end
  
  if ~isempty(xX.iC)
    fprintf('\Warning: Covariates are not recommended for repeated measures Anova and will not properly work.\n');
    fprintf('Please remove any covariates from your longitudinal design.\n\n');
  end
  
  fprintf('\nPlease note that permutation is only done within subjects for repeated Anova.\n');
else
  exch_block_labels = ones(1,n_data);
end

if ~test_mode
  % check for meshes
  if spm_mesh_detect(VY)
    mesh_detected = true;
  else
    mesh_detected = false;
  end
  
  % set E according to type data
  if job.tbss || mesh_detected
    E = 1.0; H = 2.0;
  else
    E = 0.5; H = 2.0;
  end

  % check for mask image that should exist for any analysis
  if exist(fullfile(cwd, 'mask.img'),'file')
    file_ext = '.img';
  elseif exist(fullfile(cwd, 'mask.nii'),'file')
    file_ext = '.nii';
  elseif exist(fullfile(cwd, 'mask.gii'),'file')
    file_ext = '.gii';
  else
    fprintf('WARNING: No mask file found. Switch to test mode.\n');
    test_mode = true;
  end

else
  % for test_mode
  mesh_detected = false;
end

% check whether 3D filters were selected for surfaces meshes or whether CAT12 is installed
if mesh_detected
  if ~exist('spm_cat12','file')
    error('For using surface analysis you need to install CAT12.');
  end  
end

if ~test_mode
  % get mask file
  if isempty(job.mask)
    maskname = fullfile(cwd,['mask' file_ext]);
  else
    if ~isempty(job.mask{1})
      maskname = job.mask{1};
    else    
      maskname = fullfile(cwd,['mask' file_ext]);
    end
  end
  
  % load mask
  try
    Vmask = spm_data_hdr_read(maskname);
  catch
    if mesh_detected
      maskname = spm_select(1,'mesh','select surface mask');
    else
      maskname = spm_select(1,'image','select mask image');
    end
    Vmask = spm_data_hdr_read(maskname);
  end
    
  % if first image was not found you have to select all files again
  if ~exist(VY(1).fname,'file')
  
    fprintf('Data not found. Please select data in the order defined in the design matrix.\n');
    n = size(SPM.xY.VY,1);
    if mesh_detected
      P = spm_select(n,'mesh','select surfaces');
    else
      P = spm_select(n,'image','select images');
    end
    
    VY = spm_data_hdr_read(P);
    
    %-Apply gSF to memory-mapped scalefactors to implement scaling
    %--------------------------------------------------------------------------
    for i = 1:n
      VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*SPM.xGX.gSF(i);
      if mesh_detected
        if isfield(VY(i).private.private.data{1}.data,'scl_slope')
          VY(i).private.private.data{1}.data.scl_slope = ...
              VY(i).private.private.data{1}.data.scl_slope * SPM.xGX.gSF(i);
          VY(i).private.private.data{1}.data.scl_inter = ...
              VY(i).private.private.data{1}.data.scl_inter * SPM.xGX.gSF(i);
        end
      else
        VY(i).private.dat.scl_slope = ...
            VY(i).private.dat.scl_slope * SPM.xGX.gSF(i);
        VY(i).private.dat.scl_inter = ...
            VY(i).private.dat.scl_inter * SPM.xGX.gSF(i);
      end
    end
    
    SPM.xY.VY = VY;
      
    % update SPM
    if size(SPM.xY.VY,1)==n
      save(job.data{1},'SPM','-v7.3');
    else
      fprintf('ERROR: Number of files is not correct\n');
      return
    end
  end
        
  if voxel_covariate
    % if first image of voxel-wise covariate was not found you have to select all files again
    if ~exist(SPM.xC(voxel_covariate).VC(1).fname,'file')
  
      fprintf('Covariate data not found. Please select covariate data in the order defined in the design matrix.\n');
      n = numel(SPM.xC(voxel_covariate).P);
      if mesh_detected
        P = spm_select(n,'mesh','select surfaces');
      else
        P = spm_select(n,'image','select images');
      end
      SPM.xC(voxel_covariate).P = cellstr(P);
      
      if size(P,1)==n
        save(job.data{1},'SPM','-v7.3');
      else
        fprintf('ERROR: Number of files is not correct\n');
        return
      end
    end
    
    VC = spm_data_hdr_read(char(SPM.xC(voxel_covariate).P));
    n = numel(VC);
    
    %-Apply gSF to memory-mapped scalefactors to implement scaling for
    % voxel-wise covariate
    %--------------------------------------------------------------------------
    if isfield(SPM.xC(voxel_covariate),'gSF')
      for i = 1:n
        VC(i).pinfo(1:2,:) = VC(i).pinfo(1:2,:)*SPM.xC(voxel_covariate).gSF(i);
        if mesh_detected
          if isfield(VC(i).private.private.data{1}.data,'scl_slope')
            VC(i).private.private.data{1}.data.scl_slope = ...
                VC(i).private.private.data{1}.data.scl_slope * SPM.xC(voxel_covariate).gSF(i);
            VC(i).private.private.data{1}.data.scl_inter = ...
                VC(i).private.private.data{1}.data.scl_inter * SPM.xC(voxel_covariate).gSF(i);
          end
        else
          VC(i).private.dat.scl_slope = ...
              VC(i).private.dat.scl_slope * SPM.xC(voxel_covariate).gSF(i);
          VC(i).private.dat.scl_inter = ...
              VC(i).private.dat.scl_inter * SPM.xC(voxel_covariate).gSF(i);
        end
      end
    end
    SPM.xC(voxel_covariate).VC = VC;
  end

  % check whether mask images fits to the data
  if mesh_detected, dim_index = 1; else dim_index=1:3; end
  if sum(sum((Vmask.mat-VY(1).mat).^2)) > 1e-6 || any(Vmask.dim(dim_index) ~= VY(1).dim(dim_index))
    error('Mask must have the same dimensions and orientation as the data.');
  end
  
  if voxel_covariate
    do_resample = false;
    if sum(sum((Vmask.mat-VC(1).mat).^2)) > 1e-6 || any(Vmask.dim(dim_index) ~= VC(1).dim(dim_index))
      fprintf('Covariate data has different dimensions and orientation as the other data. Thus, your covariate data will be resampled.\n');
      fprintf('Please check that your covariate data are co-registered to your other data!\n');
      do_resample = true;
    end
  end

  % read mask and data
  mask = spm_data_read(Vmask);
  
  ind_mask = find(mask>0);
  n = numel(VY);
  
  if ~isempty(ind_mask)
    Y = zeros([length(ind_mask) n],'single');
    if voxel_covariate
      C = zeros([length(ind_mask) n],'single');
      meanC = zeros(length(ind_mask),1,'single');
    else
      C = [];
    end
    
    % load data
    fprintf('Load data ')
    for i=1:n
      fprintf('.')
      tmp = spm_data_read(VY(i));
      Y(:,i) = tmp(ind_mask);
      if voxel_covariate
        
        % we need to resample voxel-wise covariate
        if do_resample
          tmp = zeros(Vmask.dim);
          for z=1:Vmask.dim(3)
            M = spm_matrix([0 0 z]);
            M1  = Vmask(1).mat\VC(i).mat\M;
            tmp(:,:,z) = spm_slice_vol(VC(i),M1, Vmask.dim(1:2),[1 0]);
          end  
        else
          tmp = spm_data_read(VC(i));
        end
        
        tmp = tmp(ind_mask);
        C(:,i) = tmp;
        meanC = meanC + tmp;
      end
    end
    fprintf('\n')
    clear tmp;
  
    % whitening matrix
    W = single(full(xX.W));
    Y = Y*W;

    % mean correct voxel-wise covariate
    if voxel_covariate
      meanC = meanC/n;
      for i=1:n
        C(:,i) = C(:,i) - meanC;
      end
      C = C*W;
    end

  else
    error('Empty mask.');
  end
  
  t0 = zeros(Vmask.dim);
  t  = zeros(Vmask.dim);
  
end % if ~test_mode

% if just one contrast is defined, use this contrast and skip interactive selection
if isfield(SPM,'xCon') && numel(SPM.xCon) == 1 && ~isfinite(Ic0)
  Ic0 = 1;
end

% interactively select contrast(s) if necessary
if numel(Ic0)==1 && ~isfinite(Ic0) && numel(SPM.xCon) > 1
  [Ic0,xCon] = spm_conman(SPM,'T&F',Inf,...
    '  Select contrast(s)...',' ',1);
  SPM.xCon = xCon;
end

% go through all contrasts
for con = 1:length(Ic0)
    
  Ic = Ic0(con);
  xCon = SPM.xCon(Ic);
  
  n_perm = job.conspec.n_perm(1);
  if numel(job.conspec.n_perm) > 1
    n_perm_break = job.conspec.n_perm(2);
  end
      
  if length(Ic) > 1
    fprintf('ERROR: No conjunction allowed.\n');
    return
  end
  
  fprintf('Use contrast #%d of %s\n',Ic,job.data{1})

  % get contrast and name
  c0 = xCon.c;  
  F_contrast_multiple_rows = 0;
  
  % for F-contrasts if rank is 1 we can use the first row
  if strcmp(xCon.STAT,'F')
    if rank(c0) == 1
      c0 = c0(:,1);
    else
      F_contrast_multiple_rows = 1;
    end
  end

  [indi, indj] = find(c0~=0);
  ind_X = unique(indi)';
  xCon.ind_X = ind_X;
  
  % check for contrasts that are defined for columns with subject effects
  if ~isempty(xX.iB)
    if max(ind_X) > min(xX.iB)
      fprintf('ERROR: No contrasts on subjects/block effects allowed.\n');
      return
    end
  end

  % find exchangeability blocks using contrasts without zero values
  exch_blocks   = c0(ind_X,:);
    
  n_exch_blocks = length(ind_X);
  
  % recognize effects of interest contrast for F-tests
  if F_contrast_multiple_rows && size(exch_blocks,2) == n_exch_blocks
    is_eoi = all(all(exch_blocks == eye(n_exch_blocks)));
    if is_eoi
      n_exch_blocks = 1;
    end
  end
  
  % check for exchangeability blocks and design matrix
  if n_exch_blocks == 1
    n_cond = length(find(xX.iH==ind_X)); % check whether the contrast is defined at columns for condition effects
  else
    n_cond = 0;
    n_data_cond = [];
    for k=1:length(xX.iH)
      n_data_cond = [n_data_cond sum(xX.X(:,xX.iH(k)))];
    end
    
    % for F-contrast with multiple rows n_cond is always n_exch_blocks
    if F_contrast_multiple_rows && length(xX.iH) > 1
      n_cond = n_exch_blocks;
    elseif F_contrast_multiple_rows && length(xX.iH) == 1
      n_cond = 0;
    else
      for j=1:n_exch_blocks
        col_exch_blocks = find(c0==exch_blocks(j));
        for k=1:length(col_exch_blocks)
          n_cond = n_cond + length(find(xX.iH==col_exch_blocks(k)));
        end
      end  
    end
    
  end

  use_half_permutations = 0;
  % check if sample size is equal for both conditions
  if n_cond == 2
    try
      % repated Anova or F-test don't allow to use only half of the permutions
      if repeated_anova || strcmp(xCon.STAT,'F')
        use_half_permutations = 0;
      elseif sum(n_data_cond(c0==exch_blocks(1))) == sum(n_data_cond(c0==exch_blocks(2)))
        use_half_permutations = 1;
      end
    end
  end
  
  ind_exch_blocks = cell(n_exch_blocks,1);
  for j=1:n_exch_blocks
    if strcmp(xCon.STAT,'T')
      ind_exch_blocks{j} = find(c0==exch_blocks(j));
    else
      ind_exch_blocks{j} = ind_X(j);
    end
  end

  fprintf('\n');
  
  % check design
  interaction_design = false;
  switch n_cond
  case 0 % correlation
    label = 1:n_data;

    % we have to correct for some F-contrasts (i.e. effects of interest
    % with eyes)
    if F_contrast_multiple_rows && is_eoi
      is_one = find(any(c0'));
      for j=1:numel(is_one)
        ind_exch_blocks{j} = is_one(j);
      end
      ind_exch_blocks = ind_exch_blocks';
    end
    
    if n_exch_blocks >= 2 && any(diff(exch_blocks(:))) % # exch_blocks >1 & differential contrast
      fprintf('Interaction design between two or more regressors found\n')
      interaction_design = true;

      % remove all entries where contrast is not defined
      % this does not work for all data CG 20200829
      % label(all(xX.X(:,ind_X)==0,2)) = [];
    else
      if repeated_anova
        fprintf('Repeated Anova with contrast for covariate found\n');
      else
        fprintf('Multiple regression design found\n');
      end
    end    
  case 1 % one-sample t-test
    fprintf('One sample t-test found\n');
    
    % use exchangeability blocks for labels
    label = zeros(1,n_data);
    for j=1:n_exch_blocks
      for k=1:length(ind_exch_blocks{j})
        label(xX.X(:,ind_exch_blocks{j}(k))~=0) = j;
      end
    end
  otherwise  % Anova with at least 2 groups
    if repeated_anova
      fprintf('Repeated Anova found\n');
    else
      fprintf('Anova found\n');
    end

    % use exchangeability blocks for labels
    label = zeros(1,n_data);
    for j=1:n_exch_blocks
      for k=1:length(ind_exch_blocks{j})
        label(xX.X(:,ind_exch_blocks{j}(k))~=0) = j;
      end
    end
  end

  fprintf('\n')

  % get index for label values > 0
  ind_label  = find(label > 0);
  
  n_data_with_contrast = length(ind_label);
  
  % estimate # of permutations
  % Anova/correlation: n_perm = (n1+n2+...+nk)!/(n1!*n2!*...*nk!)
  if n_cond ~=1  % Anova/correlation
    n_perm_full = factorial(n_data_with_contrast);
    single_subject = 0;
    
    for i=1:n_cond
      % check whether only a single subject is in one group
      if length(find(label == i)) == 1
        single_subject = 1;
      end
      n_perm_full = n_perm_full/factorial(length(find(label == i)));
    end
    
    if isnan(n_perm_full)
      % correct number of permutations for large samples when factorial is not working
      if (n_cond == 2) && (single_subject == 1)
        n_perm_full = n_data_with_contrast;
      else
        n_perm_full = realmax;
      end
    end

    % find where data are defined for that contrast
    if ~isempty(find(xX.iH == ind_X(1), 1))
      % first checking whether contrasts are defined for iH
      ind_data_defined = find(any(xX.X(:,xX.iH(ind_X)),2));
    else
      ind_data_defined = find(any(xX.X(:,ind_X),2));
    end
    
    % correct ind_label and n_data_with_contrast using ind_data_defined
    ind_label  = ind_data_defined';
    n_data_with_contrast = length(ind_label);
    
    % and restrict exchangeability block labels to those rows
    exch_block_labels_data_defined = exch_block_labels(ind_data_defined);

    % Repated Anova: n_perm = n_cond1!*n_cond2!*...*n_condk!
    % for a full model where each condition is defined for all subjects the easier
    % estimation is: n_perm = (n_cond!)^n_subj
    % check that no regression analysis inside repeated anova is used
    if repeated_anova && n_cond~=0
      n_subj = max(exch_block_labels_data_defined);
      n_perm_full = 1;
      for k=1:n_subj
        n_cond_subj = length(find(exch_block_labels_data_defined == k));
        n_perm_full = n_perm_full*factorial(n_cond_subj);
      end
    else
      n_perm_full = round(n_perm_full);
    end
    
  else  % one-sample t-test: n_perm = 2^n
    n_perm_full = 2^n_data_with_contrast;
    exch_block_labels_data_defined = exch_block_labels;
    ind_data_defined = ind_label;
  end

  % sometimes for F-tests with multiple independent rows the design cannot be fully recognized
  % and # of permutations is wrong
  if n_perm_full == 1 && F_contrast_multiple_rows
    fprintf('ERROR: This F-contrast and type of design with multiple independent rows is not yet supported.\n');
    return
  end
  
  if n_perm_full < n_perm
    fprintf('Warning: Maximum number of possible permutations is lower than defined number of permutations: %d\n',n_perm_full);
  end
  
  n_perm = min([n_perm n_perm_full]);

  fprintf('Number of permutations: %d\n',n_perm);
  
  if use_half_permutations
    fprintf('Equal sample sizes: Use half the number of permutations.\n');
  end
  
  fprintf('Exchangeability block/variable: ');
  fprintf('%d ',unique(cell2mat(ind_exch_blocks)));
  fprintf('\n');
  fprintf('# of conditions: %d\n',n_cond);
           
  % Guttman partioning of design matrix into effects of interest X and nuisance variables Z
  X = xX.X(:,ind_X);
  ind_Z = [xX.iH xX.iC xX.iB xX.iG];
  ind_Z(ind_X) = [];
  Z = xX.X(:,ind_Z);
    
  Hz = Z*pinv(Z);
  Rz = eye(size(X,1)) - Hz;

  % if Hz is zero or Ic is empty then no confounds were found and we can skip the time-consuming
  % Freedman-Lane permutation
  if (all(~any(Hz)) || isempty(xX.iC)) || all(~any(diff(Hz))) || (interaction_design && numel(xX.iC) == numel(ind_X))
    exist_nuisance = false;
  else
    exist_nuisance = true;
  end
  
  if ~exist_nuisance && nuisance_method > 0
    fprintf('No nuisance variables were found: Use Draper-Stoneman permutation.\n\n');
    nuisance_method = 0;
  end

  if nuisance_method > 0 && repeated_anova
    fprintf('Use Draper-Stoneman permutation for repeated measures Anova.\n\n');
    nuisance_method = 0;
  end

  switch nuisance_method 
  case 0
    str_permutation_method = 'Draper-Stoneman';
  case 1
    str_permutation_method = 'Freedman-Lane';
  case 2
    str_permutation_method = 'Smith';
  end

  % name of contrast
  c_name0 = deblank(xCon.name);

  if test_mode
    c_name = '';
  else
    c_name = sprintf('%s (E=%1.1f H=%1.1f %s) ',c_name0, E, H, str_permutation_method);
  end

  if ~test_mode
            
    % compute unpermuted t/F-map
    if voxel_covariate

      % check which pinv-method is faster and use that one for all permutations   
      X = xX.W*xX.X;
  
      tstart1 = tic;
      for i=1:20, pX = pinv(X); end
      telapsed1 = toc(tstart1);

      tstart2 = tic;
      for i=1:20, pX = pinv2(X); end
      telapsed2 = toc(tstart2);

      if 1.1*telapsed2 < telapsed1
        pinv_method = 2;
        fprintf('Use faster pinv2 function\n');
      else
        pinv_method = 1;
      end
      
      [t0, df2] = calc_GLM_voxelwise(Y,xX,SPM.xC(voxel_covariate),xCon,ind_mask,VY(1).dim,C,[],ind_X,pinv_method);
    else
      [t0, df2] = calc_GLM(Y,xX,xCon,ind_mask,VY(1).dim);
    end

    df1 = size(xCon.c,2);
    
    % transform to z statistic
    if convert_to_z
      % use faster z-transformation of SPM for T-statistics
      if strcmp(xCon.STAT,'T')
        t0 = spm_t2z(t0,df2);
      else
        t0 = palm_gtoz(t0,df1,df2);
      end
    end
    
    mask_0   = (t0 == 0);
    mask_1   = (t0 ~= 0);
    mask_P   = (t0 > 0);
    mask_N   = (t0 < 0);
    mask_NaN = (mask == 0);
    found_P  = sum(mask_P(:)) > 0;
    found_N  = sum(mask_N(:)) > 0;
      
    % remove all NaN and Inf's
    t0(isinf(t0) | isnan(t0)) = 0;
  
    % sometimes z-transformation produces neg. values even for F-statistics
    if strcmp(xCon.STAT,'F')
      t0(t0 < 0) = 0;
    end
  
    % get parametric p-values for comparison
    tname = sprintf('spm%s_%04d',xCon.STAT,Ic);
    tname = fullfile(cwd,[tname file_ext]);

    if ~exist(tname,'file') && ~voxel_covariate
      spm_contrasts(SPM,Ic);
    end

    if ~voxel_covariate
      Z0 = spm_data_read(tname);

      Pt = zeros(size(Z0));
      if strcmp(xCon.STAT,'T')
        if found_P
          Pt(mask_P) = 1-spm_Tcdf(Z0(mask_P),df2);
        else
          Pt(mask_N) = spm_Tcdf(Z0(mask_N),df2)-1;
        end
      else
        if found_P
          Pt(mask_P) = 1-spm_Fcdf(Z0(mask_P),[df1, df2]);
        else
          Pt(mask_N) = spm_Fcdf(Z0(mask_N),[df1, df2])-1;
        end
      end
          
      % Check correlation between parametric and non-parametric T/F-values.
      % Low correlation points to issues with image mask and we have to use a
      cc = corrcoef(Z0(:),t0(:));
      mask_shared = Z0 ~= 0;
      if cc(1,2) < 0.85 && isempty(job.mask)
        % check whether mask size differes and create a shared mask 
        if sum(Z0(:) ~= 0) ~= sum(t0(:) ~= 0)
          fprintf('\nWARNING: Large discrepancy between parametric and non-parametric statistic found (cc=%g) which either points to different image masks or to missing absolute threshold for VBM analysis.\n',cc(1,2));
          mask_shared = Z0 ~= 0 & t0 ~=0;
        else
          fprintf('\nWARNING: Large discrepancy between parametric and non-parametric statistic found (cc=%g) which is likely due to creating parametric statistics in fMRI mode, which slightly handles noise differently.\n',cc(1,2));
        end
      end
          
    else
      mask_shared = ones(size(t0));
    end

    % get dh for unpermuted map
    dh = max(abs(t0(:)))/n_steps_tfce;
  
    % calculate tfce of unpermuted t-map
    if mesh_detected
      if ~isa(SPM.xVol.G,'gifti')
        % check whether path is correct and file exist
        if ~exist(SPM.xVol.G,'file')
          [pathG,nameG,extG] = spm_fileparts(SPM.xVol.G);
          % use new path
          if ~isempty(strfind(pathG,'_32k'))
            SPM.xVol.G = fullfile(fileparts(which('cat12')),'templates_surfaces_32k',[nameG extG]);
          else
            SPM.xVol.G = fullfile(fileparts(which('cat12')),'templates_surfaces',[nameG extG]);
          end
        end
        SPM.xVol.G = gifti(SPM.xVol.G);
      end
      tfce0 = tfce_mesh(SPM.xVol.G.faces, t0, dh, E, H)*dh;
    else
      
      % measure computation time to test whether multi-threading causes issues
      % start with single-threading for unpermuted data
      tstart = tic;
      % only estimate neg. tfce values for non-positive t-values
      if found_N
        tfce0 = tfceMex_pthread(t0,dh,E,H,1,1)*dh;
      else
        tfce0 = tfceMex_pthread(t0,dh,E,H,0,1)*dh;
      end
      telapsed = toc(tstart);
    end

    % prepare output files
    Vt = VY(1);
    Vt.dt(1)    = 16;
    Vt.pinfo(1) = 1;

    %---------------------------------------------------------------
    % save unpermuted t map
    %---------------------------------------------------------------
    name = sprintf('%s_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d %s',xCon.STAT,Ic, str_permutation_method);
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,t0);
  
    %---------------------------------------------------------------
    % save unpermuted TFCE map
    %---------------------------------------------------------------
    name = sprintf('TFCE_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d %s',Ic, str_permutation_method);
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,tfce0);
  
    % get largest tfce
    tfce0_max = max(tfce0(:));
    t0_max    = max(t0(:));
    tfce0_min = min(tfce0(:));
    t0_min    = min(t0(:));
        
    % prepare countings
    tperm        = zeros(size(t));
    tfceperm     = zeros(size(t));
    t_min        = [];
    t_max        = [];
    t_max_th     = [];
    tfce_min     = [];
    tfce_max     = [];
    tfce_max_th  = [];
  
  end % test_mode

  % general initialization
  try % use try commands to allow batch mode without graphical output
    Fgraph = spm_figure('GetWin','Graphics');
    spm_figure('Clear',Fgraph);
    figure(Fgraph)
  
    h = axes('position',[0.45 0.95 0.1 0.05],'Units','normalized','Parent',...
      Fgraph,'Visible','off');
      
    text(0.5,0.6,c_name,...
      'FontSize',spm('FontSize',10),...
      'FontWeight','Bold',...
      'HorizontalAlignment','Center',...
      'VerticalAlignment','middle')
  
    text(0.5,0.25,spm_str_manip(spm_fileparts(job.data{1}),'a80'),...
      'FontSize',spm('FontSize',8),...
      'HorizontalAlignment','Center',...
      'VerticalAlignment','middle')
  end
  
  % check that label has correct dimension
  sz = size(label);
  if sz(1)>sz(2)
    label = label';
  end
    
  stopStatus = false;
  if ~test_mode, tfce_progress('Init',n_perm,'Calculating','Permutations'); end
  
  % update interval for progress bar
  progress_step = max([1 round(n_perm/100)]);

  % Regression design found where contrast is defined for covariate?
  if ~isempty(xX.iC) && all(ismember(ind_X,SPM.xX.iC))
    ind_label_gt0 = find(label(ind_data_defined) > 0);
  else
    ind_label_gt0 = find(label > 0);
  end
  
  unique_labels = unique(label(ind_label_gt0));
  n_unique_labels = length(unique_labels);
  
  perm = 1;
  check_validity = false;
  while perm<=n_perm

    % randomize subject vector
    if perm==1 % first permutation is always unpermuted model
      if n_cond == 1 % one-sample t-test
        rand_label = ones(1,n_data_with_contrast);
        label_matrix = rand_label;
      else % correlation or Anova
        rand_order = ind_label;
        rand_order_sorted = rand_order;
        label_matrix = rand_order;
      end
    else
      % init permutation and
      % check that each permutation is used only once
      if n_cond == 1 % one-sample t-test
        rand_label = sign(randn(1,n_data_with_contrast));
        while any(ismember(label_matrix,rand_label,'rows'))
          rand_label = sign(randn(1,n_data_with_contrast));
        end
      else % correlation or Anova
        
        % permute inside exchangeability blocks only
        rand_order = zeros(1,n_data_with_contrast);
        rand_order_sorted = zeros(1,n_data_with_contrast);
        for k = 1:max(exch_block_labels_data_defined)
          ind_block   = find(exch_block_labels_data_defined == k);
          n_per_block = length(ind_block);
          rand_order(ind_block) = ind_label(ind_block(randperm(n_per_block)));
        end
        
        % go through defined labels and sort inside
        for k=1:n_unique_labels
          ind_block = find(label(ind_label_gt0) == unique_labels(k));
          rand_order_sorted(ind_block) = sort(rand_order(ind_block));
        end

        % check whether this permutation was already used
        count_trials = 0;
        while any(ismember(label_matrix,rand_order_sorted,'rows'))
          count_trials = count_trials + 1;
          
          % stop the permutation loop for too many successless trials for finding 
          % new permutations
          if count_trials > 100000
            fprintf('Stopped after %d permutations because there were too many successless trials for finding new permutations.\n',perm);
            fprintf('Probably there are some missing values for some subjects and the number of maximal permutations was too high.\n');
            n_perm = perm; % stop the permutation loop
            break
          end
          
          for k = 1:max(exch_block_labels_data_defined)
            ind_block   = find(exch_block_labels_data_defined == k);
            n_per_block = length(ind_block);
            rand_order(ind_block) = ind_label(ind_block(randperm(n_per_block)));
          end
          
          % go through defined labels and sort inside
          for k=1:n_unique_labels
            ind_block = find(label(ind_label_gt0) == unique_labels(k));
            rand_order_sorted(ind_block) = sort(rand_order(ind_block));
          end

        end
      end    
    end   
    
    % create permutation set
    Pset = sparse(n_data,n_data);
    if n_cond == 1 % one-sample t-test
      for k=1:n_data_with_contrast
        Pset(ind_label(k),ind_label(k)) = rand_label(k);
      end
    else % correlation or Anova
      for k=1:n_data_with_contrast
        Pset(rand_order_sorted(k),ind_label(k)) = 1;
      end
    end

    % add Stop button after 20 iterations
    try % use try commands to allow batch mode without graphical output
      if perm==21
        hStopButton = uicontrol(Fgraph,...
          'position',[10 10 70 20],...
          'style','toggle',...
          'string','Stop',...
          'backgroundcolor',[1 .5 .5]); % light-red
      end
    
      if perm>=21
        stopStatus = get(hStopButton,'value');
      end
    
      % check Stop status
      if (stopStatus == true)
        fprintf('Stopped after %d iterations.\n',perm);
        break; % stop the permutation loop
      end
    end
      
    % change design matrix according to permutation order
    % only permute columns, where contrast is defined
    Xperm = xX.X;
    Xperm_debug = xX.X;
    Wperm = xX.W;

    switch nuisance_method 
    case 0 % Draper-Stoneman is permuting X
      Xperm(:,ind_X) = Pset*Xperm(:,ind_X);
%      if n_cond ~= 1
%        Wtmp = Pset*xX.W;
%        Wperm(ind_data_defined,ind_data_defined) = Wtmp(ind_data_defined,ind_data_defined);
%      end
    case 1 % Freedman-Lane is permuting Y
      Xperm = xX.X;
    case 2 % Smith method is additionally orthogonalizing X with respect to Z
      Xperm(:,ind_X) = Pset*Rz*Xperm(:,ind_X);
%      if n_cond ~= 1
%        Wtmp = Pset*Rz*xX.W;
%        Wperm(ind_data_defined,ind_data_defined) = Wtmp(ind_data_defined,ind_data_defined);
%      end
    end
            
    Xperm_debug(:,ind_X) = Pset*Xperm_debug(:,ind_X);
    
    % correct interaction designs
    % # exch_blocks >1 & # cond == 0 & differential contrast
    if n_exch_blocks >= 2 && n_cond==0 && ~all(exch_blocks(:))
      Xperm2 = Xperm;
      Xperm2(:,ind_X) = 0;
      for j=1:n_exch_blocks
        ind_Xj = find(xX.X(:,ind_X(j)));
        Xperm2(ind_Xj,ind_X(j)) = sum(Xperm(ind_Xj,ind_X),2);
      end
      Xperm = Xperm2;

      Xperm_debug2 = Xperm_debug;
      Xperm_debug2(:,ind_X) = 0;
      for j=1:n_exch_blocks
        ind_Xj = find(xX.X(:,ind_X(j)));
        Xperm_debug2(ind_Xj,ind_X(j)) = sum(Xperm_debug(ind_Xj,ind_X),2);
      end
      Xperm_debug = Xperm_debug2;
    end
    
    if show_permuted_designmatrix
      % scale covariates and nuisance variables to a range 0.8..1
      % to properly display these variables with indicated colors
      if ~isempty(xX.iC)
        val = Xperm_debug(:,xX.iC);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        Xperm_debug(:,xX.iC) = val;
      end
      if ~isempty(xX.iG)
        val = Xperm_debug(:,xX.iG);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        Xperm_debug(:,xX.iG) = val;
      end
      if ~isempty(xX.iH) && n_cond==1 % one-sample t-test
        val = Xperm_debug(:,xX.iH);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        Xperm_debug(:,xX.iH) = val;
      end
      
      % use different colors for indicated columns
      Xperm_debug(:,xX.iH) = 16*Xperm_debug(:,xX.iH);
      Xperm_debug(:,xX.iC) = 24*Xperm_debug(:,xX.iC);
      Xperm_debug(:,xX.iB) = 32*Xperm_debug(:,xX.iB);
      Xperm_debug(:,xX.iG) = 48*Xperm_debug(:,xX.iG);

      if n_cond==1 % one-sample t-test
        for j=1:n_data_with_contrast
          if rand_label(j) > 0
            Xperm_debug(ind_label(j),ind_X) = 60*rand_label(j)*Xperm_debug(ind_label(j),ind_X);
          else
            Xperm_debug(ind_label(j),ind_X) = 56*rand_label(j)*Xperm_debug(ind_label(j),ind_X);
          end
        end
      else % correlation or Anova
        % scale exchangeability blocks also to values 0.8..1
        val = Xperm_debug(:,ind_X);
        ind0 = (val==0);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
      
        % rescue zero entries
        val(ind0) = 0;
      
        Xperm_debug(:,ind_X) = 60*val;
      end

    end
          
    show_plot = 0;
    if use_half_permutations
      if ~rem(perm,progress_step) || ~rem(perm+1,progress_step)
        show_plot = 1;
      end
    else
      if ~rem(perm,progress_step)
        show_plot = 1;
      end
    end

    % display permuted design matrix
    try
      if show_permuted_designmatrix && show_plot
        figure(Fgraph);
        subplot(2,2,3);
        image(Xperm_debug); axis off
        title('Permuted design matrix','FontWeight','bold');
      
        % use different colormap for permuted design matrix
        cmap = jet(64);
      
        % zero values should be always black
        cmap(1,:) = [0 0 0];
        colormap(cmap)
      
        % show legend only once
        if perm <= progress_step
          subplot(2,2,4); axis off
        
          % color-coded legend
          y = 1.0;
          text(-0.2,y, 'Columns of design matrix: ', 'Color',cmap(1, :),'FontWeight','Bold','FontSize',10); y = y - 0.10;
          text(-0.2,y,['Exch. block: ' num2str_short(unique(cell2mat(ind_exch_blocks))')], 'Color',cmap(60,:),'FontWeight','Bold','FontSize',10); y = y - 0.05;
          if ~isempty(xX.iH)
            text(-0.2,y, ['iH - Indicator variable: ' num2str_short(xX.iH)], 'Color',cmap(16,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05; 
          end
          if ~isempty(xX.iC)
            text(-0.2,y, ['iC - Covariate: ' num2str_short(xX.iC)], 'Color',cmap(24,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
          if ~isempty(xX.iB)
            text(-0.2,y, ['iB - Block variable: ' num2str_short(xX.iB)], 'Color',cmap(32,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
          if ~isempty(xX.iG)
            text(-0.2,y, ['iG - Nuisance variable: ' num2str_short(xX.iG)], 'Color',cmap(48,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
        end
      end
    end
    
    if ~test_mode
      % calculate permuted t-map
      if perm == 1
        t    = t0;
        tfce = tfce0;
        % prepare null distribution
        if save_null_distribution
          null_distribution = zeros(size(t));
        end
      else
        xXperm   = xX;
        xXperm.X = Xperm;        
        xXperm.W = Wperm;

        % Freedman-Lane permutation of data
        if nuisance_method == 1
          t = calc_GLM(Y*(Pset'*Rz),xXperm,xCon,ind_mask,VY(1).dim);
        else
          if voxel_covariate
            t = calc_GLM_voxelwise(Y,xXperm,SPM.xC(voxel_covariate),xCon,ind_mask,VY(1).dim,C,Pset,ind_X,pinv_method);
          else
            t = calc_GLM(Y,xXperm,xCon,ind_mask,VY(1).dim);
          end
        end

        if convert_to_z
          % use faster z-transformation of SPM for T-statistics
          if strcmp(xCon.STAT,'T')
            t(mask_1) = spm_t2z(t(mask_1),df2);
          else
            t(mask_1) = palm_gtoz(t(mask_1),df1,df2);
          end
        end

        % update null-distribution
        if save_null_distribution
          null_distribution(mask_1) = null_distribution(mask_1) + t(mask_1);
        end
        
        % remove all NaN and Inf's
        t(isinf(t) | isnan(t)) = 0;
        
        % use individual dh
        dh = max(abs(t(:)))/n_steps_tfce;
        
        % compute tfce
        if mesh_detected
          tfce = tfce_mesh(SPM.xVol.G.faces, t, dh, E, H)*dh;
        else
          
          % measure computation time for 1st permutation to test whether multi-threading causes issues
          if perm==3 && ~singlethreaded, tstart = tic; end
          
          % only estimate neg. tfce values for non-positive t-values
          if min(t(:)) < 0
            tfce = tfceMex_pthread(t,dh,E,H,1,singlethreaded)*dh;
          else
            tfce = tfceMex_pthread(t,dh,E,H,0,singlethreaded)*dh;
          end
          
          % if multi-threading takes 3x longer then force single-threading
          % because for some unknown reason multi-threading is not working properly
          if perm==3 && ~singlethreaded
            telapsed_multi = toc(tstart);
            if (telapsed_multi > 3*telapsed)
              fprintf('Warning: Multi-threading disabled because of run-time issues.\n');
              singlethreaded = 1;
            end
          end

        end
        
      end
      
      % use (too liberal) method for estimating maximum statistic from old release 
      % r184 for compatibility purposes only that was estimating max/min statistics
      % only inside pos./neg. effects and not both
      if old_method_stat
        mask_stat_P = mask_P;
        mask_stat_N = mask_N;
      else  
        mask_stat_P = mask_1;
        mask_stat_N = mask_1;
      end
    end % test_mode
    
    % update label_matrix to check for unique permutations
    if use_half_permutations
      if perm>1
        label_matrix = [label_matrix; rand_order_sorted; [rand_order_sorted(label(ind_label) == 2) rand_order_sorted(label(ind_label) == 1)]];
      end
      if ~test_mode
        % maximum statistic
        t_max    = [t_max    max(t(mask_stat_P))    -min(t(mask_stat_N))];
        t_min    = [t_min    min(t(mask_stat_N))    -max(t(mask_stat_P))];
        tfce_max = [tfce_max max(tfce(mask_stat_P)) -min(tfce(mask_stat_N))];
        tfce_min = [tfce_min min(tfce(mask_stat_N)) -max(tfce(mask_stat_P))];
        tperm(mask_P)    = tperm(mask_P) + 2*(t(mask_P) >= t0(mask_P));
        tperm(mask_N)    = tperm(mask_N) - 2*(t(mask_N) <= t0(mask_N));
        tfceperm(mask_P) = tfceperm(mask_P) + 2*(tfce(mask_P) >= tfce0(mask_P));
        tfceperm(mask_N) = tfceperm(mask_N) - 2*(tfce(mask_N) <= tfce0(mask_N));
        
      end
    else
      if perm>1
        if n_cond == 1 % one-sample t-test
          label_matrix = [label_matrix; rand_label];
        else
          label_matrix = [label_matrix; rand_order_sorted];
        end
      end
      if ~test_mode
        % maximum statistic
        t_max    = [t_max    max(t(mask_stat_P))];
        t_min    = [t_min    min(t(mask_stat_N))];
        tfce_max = [tfce_max max(tfce(mask_stat_P))];
        tfce_min = [tfce_min min(tfce(mask_stat_N))];
        tperm(mask_P)    = tperm(mask_P) + (t(mask_P) >= t0(mask_P));
        tperm(mask_N)    = tperm(mask_N) - (t(mask_N) <= t0(mask_N));
        tfceperm(mask_P) = tfceperm(mask_P) + (tfce(mask_P) >= tfce0(mask_P));
        tfceperm(mask_N) = tfceperm(mask_N) - (tfce(mask_N) <= tfce0(mask_N));
  
      end
    end
      
    if ~test_mode
      % use cummulated sum to find threshold
      stfce_max = sort(tfce_max);
      st_max    = sort(t_max);
  
      % find corrected thresholds
      ind_max  = ceil((1-alpha).*length(st_max));
      t_max_th = [t_max_th; st_max(ind_max)];
      if use_half_permutations
        t_max_th = [t_max_th; st_max(ind_max)];
      end
      
      ind_max     = ceil((1-alpha).*length(stfce_max));
      tfce_max_th = [tfce_max_th; stfce_max(ind_max)];
      if use_half_permutations
        tfce_max_th = [tfce_max_th; stfce_max(ind_max)];
      end
      
      % plot thresholds and histograms      
      try
        if show_plot
          figure(Fgraph);
          axes('position',[0 0 1 0.95],'Parent',Fgraph,'Visible','off');
          plot_distribution(stfce_max, tfce_max_th, 'tfce', alpha, col, 1, tfce0_max, tfce0_min);
          if ~show_permuted_designmatrix
            plot_distribution(st_max, t_max_th, 't-value', alpha, col, 2, t0_max, t0_min);
          end
        end
      end
    
      if numel(job.conspec.n_perm) > 1
        if perm > n_perm_break
          if isempty(find(tfce0_max > tfce_max_th(50:end,1), 1))
            fprintf('No FWE-corrected suprathreshold value after %d permutations found\n', n_perm_break);
            perm = n_perm;
          end
        end  
      end
      
      save_results = 1;
      % wait until 50 permutations are finished and skip that if voxel-wise covariate is used 
      if ~voxel_covariate && (perm > 50) 
        % after defined number of permutations check whether maximum value exceed 95% of threshold
        if (perm >= stop_if_no_FWEeffects_found) && (tfce0_max < 0.95*stfce_max(ind_max(alpha==0.05)) && -tfce0_min < 0.95*stfce_max(ind_max(alpha==0.05)))
          fprintf('Stop estimation because after %d permutations because threshold could not be exceeded.\n',perm);
          save_results = 0;
          break; % stop the permutation loop
        end
      end

      % after 500 permutations or at n_perm compare uncorrected p-values with permutations with parametric
      % p-values to check wheter something went wrong    
      % use odd numbers to consider parameter use_half_permutations
      % skip that check for voxel-wise covariates
      
      if ~voxel_covariate && ((perm == 501) || (perm >= n_perm-1)) && ~check_validity && (found_P || found_N)
        
        % estimate p-values
        nPt = tperm/perm;
    
        % check correlation between parametric and non-parametric p-values
        % exclude Pt==0.5 and Pt==1 values that can distort masked analysis values
        if found_P
          cc = corrcoef(nPt(mask_P & Pt ~=0.5 & Pt ~=1 & mask_shared),Pt(mask_P & Pt ~=0.5 & Pt ~=1 & mask_shared));
        else
          cc = corrcoef(nPt(mask_N & Pt ~=0.5 & Pt ~=1 & mask_shared),Pt(mask_N & Pt ~=0.5 & Pt ~=1 & mask_shared));
        end

        % check for low correlation between non-parametric and permutation test
        % skip check for voxel-wise covariate
        if cc(1,2) < 0.85
          % check correlation between parametric and non-parametric statistic ofr Smith or Freedman-Lane correction
          if nuisance_method > 0 
            spm('alert!',sprintf('WARNING: Large discrepancy between parametric and non-parametric statistic found! Please try a different method to deal with nuisance parameters.\n'),'',spm('CmdLine'),0);
            fprintf('\nWARNING: Large discrepancy between parametric and non-parametric statistic found (cc=%g)! Please try a different method to deal with nuisance parameters.\n',cc(1,2));
          else
            spm('alert!',sprintf('WARNING: Large discrepancy between parametric and non-parametric statistic found! Probably your design was not correctly recognized.\n'),'',spm('CmdLine'),0);
            fprintf('\nWARNING: Large discrepancy between parametric and non-parametric statistic found (cc=%g)! Probably your design was not correctly recognized.\n',cc(1,2));
          end
        else
          fprintf('\nCorrelation between between parametric and non-parametric statistic is cc=%g, which means that your design and optionally your nuisance parameters were correctly recognized.\n',cc(1,2));
        end
        check_validity = true;
      end

    end % test_mode

    
    if show_plot
      if ~test_mode, tfce_progress('Set',perm,Fgraph); end
      drawnow
    end
      
    if use_half_permutations  
      perm = perm + 2;
    else
      perm = perm + 1;
    end
  
  end
  
  if ~test_mode, tfce_progress('Clear',Fgraph); end
  
  try
    delete(hStopButton);
    spm_print;
  end
  
  if ~test_mode && save_results
    % get correct number of permutations in case that process was stopped
    n_perm = length(tfce_max);
  
    %---------------------------------------------------------------
    % corrected threshold based on permutation distribution
    %---------------------------------------------------------------
  
    % prepare output files
    Vt = VY(1);
    Vt.dt(1)    = 16;
    Vt.pinfo(1) = 1;
  
    % save ascii file with number of permutations
    name = sprintf('%s_%04d',xCon.STAT,Ic);
    fid = fopen(fullfile(cwd,[name '.txt']),'w');
    fprintf(fid,'%d\n',n_perm);
    fclose(fid);
  
    %---------------------------------------------------------------
    % save uncorrected p-values for TFCE
    %---------------------------------------------------------------
    fprintf('Save uncorrected p-values.\n');
  
    name = sprintf('TFCE_log_p_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d %s',Ic, str_permutation_method);
      
    % estimate p-values
    nPtfce = tfceperm/n_perm;
    nPtfcelog10 = zeros(size(tfce0));
  
    if found_P
      nPtfcelog10(mask_P) = -log10(nPtfce(mask_P));
    end
    
    if found_N
      nPtfce(mask_N) = -nPtfce(mask_N);
      nPtfcelog10(mask_N) =  log10(nPtfce(mask_N));
    end
    
    nPtfce(mask_0)   = 1;
    nPtfce(mask_NaN) = NaN;
    nPtfcelog10(mask_0)   = 0;
    nPtfcelog10(mask_NaN) = NaN;
  
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,nPtfcelog10);
  
    %---------------------------------------------------------------
    % save uncorrected p-values for T
    %---------------------------------------------------------------
    name = sprintf('%s_log_p_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d %s',xCon.STAT,Ic, str_permutation_method);
  
    nPtlog10 = zeros(size(t0));
  
    % estimate p-values
    nPt = tperm/n_perm;
   
    if found_P
      nPtlog10(mask_P) = -log10(nPt(mask_P));
    end
    
    if found_N
      nPt(mask_N) = -nPt(mask_N);
      nPtlog10(mask_N) =  log10(nPt(mask_N));
    end
    
    nPt(mask_0)   = 1;
    nPt(mask_NaN) = NaN;
    nPtlog10(mask_0)   = 0;
    nPtlog10(mask_NaN) = NaN;
  
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,nPtlog10);
  
    %---------------------------------------------------------------
    % save corrected p-values for TFCE
    %---------------------------------------------------------------
    fprintf('Save corrected p-values.\n');

    name = sprintf('TFCE_log_pFWE_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d FWE %s',Ic, str_permutation_method);
  
    corrP = ones(size(t));
  
    if use_gamma_tail_approximation
      fprintf('Using tail approximation from the Gamma distribution for corrected P-values.\n');
      
      if tail_approximation_wo_unpermuted_data
        ind_tail = 2:n_perm;
      else
        ind_tail = 1:n_perm;
      end
      
      if found_P
        [mu,s2,gamm1] = palm_moments(tfce_max(ind_tail)');
        corrP(mask_P) = palm_gamma(tfce0(mask_P),mu,s2,gamm1,false,1/n_perm);
      end
      
      if found_N
        [mu,s2,gamm1] = palm_moments(-tfce_min(ind_tail)');
        corrP(mask_N) = -palm_gamma(-tfce0(mask_N),mu,s2,gamm1,false,1/n_perm);
      end
      
    else

      if found_P
        for t2 = tfce_max
          %-FWE-corrected p is proportion of randomisation greater or
          % equal to statistic.
          %-Use a > b -tol rather than a >= b to avoid comparing
          % two reals for equality.
          corrP(mask_P) = corrP(mask_P) + (t2 > (tfce0(mask_P)  - tol));
        end
      end
      
      if found_N
        for t2 = tfce_min
          %-FWE-corrected p is proportion of randomisation greater or
          % equal to statistic.
          %-Use a > b -tol rather than a >= b to avoid comparing
          % two reals for equality.
          corrP(mask_N) = corrP(mask_N) - (t2 < (tfce0(mask_N) + tol));
        end
      end
      
      corrP = corrP/n_perm;  
    end
    corrPlog10 = zeros(size(tfce0));
  
    if found_P
      corrPlog10(mask_P) = -log10(corrP(mask_P));
    end
    
    if found_N
      corrP(mask_N) = -corrP(mask_N);
      corrPlog10(mask_N) =  log10(corrP(mask_N));
    end
    
    corrP(mask_0)   = 1;
    corrP(mask_NaN) = NaN;
    corrPlog10(mask_0)   = 0;
    corrPlog10(mask_NaN) = NaN;
  
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,corrPlog10);
  
    %---------------------------------------------------------------
    % save corrected p-values for T
    %---------------------------------------------------------------
    name = sprintf('%s_log_pFWE_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d FWE %s',xCon.STAT,Ic, str_permutation_method);
  
    corrP = zeros(size(t));
  
    if use_gamma_tail_approximation
      
      if found_P
        [mu,s2,gamm1] = palm_moments(t_max(ind_tail)');
        corrP(mask_P) = palm_gamma(t0(mask_P),mu,s2,gamm1,false,1/n_perm);
      end
      
      if found_N
        [mu,s2,gamm1] = palm_moments(-t_min(ind_tail)');
        corrP(mask_N) = -palm_gamma(-t0(mask_N),mu,s2,gamm1,false,1/n_perm);
      end
      
    else
      if found_P
        for t2 = t_max
          %-FWE-corrected p is proportion of randomisation greater or
          % equal to statistic.
          %-Use a > b -tol rather than a >= b to avoid comparing
          % two reals for equality.
          corrP(mask_P) = corrP(mask_P) + (t2 > t0(mask_P)  - tol);
        end
      end
      
      if found_N
        for t2 = t_min
          %-FWE-corrected p is proportion of randomisation greater or
          % equal to statistic.
          %-Use a > b -tol rather than a >= b to avoid comparing
          % two reals for equality.
          corrP(mask_N) = corrP(mask_N) - (t2 < t0(mask_N) + tol);
        end
      end
      
      corrP = corrP/n_perm;  
    end

    corrPlog10 = zeros(size(tfce0));
  
    if found_P
      corrP(mask_P) = corrP(mask_P);
      corrPlog10(mask_P) = -log10(corrP(mask_P));
    end
    
    if found_N
      corrP(mask_N) = -corrP(mask_N);
      corrPlog10(mask_N) =  log10(corrP(mask_N));
    end
  
    corrP(mask_0)   = 1;
    corrP(mask_NaN) = NaN;
    corrPlog10(mask_0)   = 0;
    corrPlog10(mask_NaN) = NaN;
  
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,corrPlog10);
  
    %---------------------------------------------------------------
    % save corrected FDR-values for TFCE
    %---------------------------------------------------------------
    fprintf('Save corrected FDR-values.\n');
  
    name = sprintf('TFCE_log_pFDR_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d FDR %s',Ic, str_permutation_method);
  
    corrPfdr = NaN(size(t));
    corrPfdrlog10 = zeros(size(tfce0));
  
    if found_P
      [snP_pos,I_pos] = sort(nPtfce(mask_P));
      if ~isempty(snP_pos)
        corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
        corrPfdr_pos(I_pos) = corrPfdr_pos;
        corrPfdr(mask_P) = corrPfdr_pos;
        corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
      end
    end
  
    if found_N
      [snP_neg,I_neg] = sort(nPtfce(mask_N));
      if ~isempty(snP_neg)
        corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
        corrPfdr_neg(I_neg) = corrPfdr_neg;
        corrPfdr(mask_N) = corrPfdr_neg;
        corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
      end
    end
      
    corrPfdrlog10(mask_0)   = 0;
    corrPfdrlog10(mask_NaN) = NaN;
  
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,corrPfdrlog10);
  
    %---------------------------------------------------------------
    % save corrected FDR-values for T
    %---------------------------------------------------------------
    name = sprintf('%s_log_pFDR_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d FDR %s',xCon.STAT,Ic, str_permutation_method);
  
    corrPfdr = NaN(size(t));
    corrPfdrlog10 = zeros(size(tfce0));
  
    if found_P
      if ~isempty(snP_pos)
        [snP_pos,I_pos] = sort(nPt(mask_P));
        corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
        corrPfdr_pos(I_pos) = corrPfdr_pos;
        corrPfdr(mask_P) = corrPfdr_pos;
        corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
      end
    end
  
    if found_N
      [snP_neg,I_neg] = sort(nPt(mask_N));
      if ~isempty(snP_neg)
        corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
        corrPfdr_neg(I_neg) = corrPfdr_neg;
        corrPfdr(mask_N) = corrPfdr_neg;
        corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
      end
    end
  
    corrPfdrlog10(mask_0)   = 0;
    corrPfdrlog10(mask_NaN) = NaN;
  
    Vt = spm_data_hdr_write(Vt);
    spm_data_write(Vt,corrPfdrlog10);

    %---------------------------------------------------------------
    % save null distribution
    %---------------------------------------------------------------
    if save_null_distribution
      fprintf('Save null distribution.\n');
      null_distribution = null_distribution/sqrt(n_perm - 1);
      name = sprintf('null%s_%04d',xCon.STAT,Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      Vt.descrip = sprintf('Null%s %04d %s',xCon.STAT,Ic, str_permutation_method);
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,null_distribution);
    end
      
  end % test_mode
    
end

colormap(gray)

%---------------------------------------------------------------
function plot_distribution(val_max,val_th,name,alpha,col,order,val0_max,val0_min)

corr = 1;

n = length(val_th);
sz_val_max = length(val_max);

% allow other thresholds depending on # of permutations
n_alpha = 3;
if sz_val_max < 1000, n_alpha = 2; end
if sz_val_max <  100, n_alpha = 1; end

% with 20 values we have the lowest possible alpha of 0.05
if sz_val_max >= 20
  alpha  = alpha(1:n_alpha);
  val_th = val_th(:,1:n_alpha);

  [hmax, xmax] = hist(val_max, 100);
      
  subplot(2,2,(2*order)-1)
  
  h = bar(xmax,hmax);
  set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);

  max_h = max(hmax);
  lim_x = xlim;

  % plot maximum observed value for unpermuted model
  hl = line([val0_max val0_max], [0 max_h]);
  set(hl,'Color',[0.3333 1 0],'LineWidth',2);
  text(0.95*val0_max,0.95*max_h,'Max. observed value ',...
    'Color',[0.3333 1 0],'HorizontalAlignment','Right','FontSize',8)
    
  % plot sign-flipped minimum observed value for unpermuted model
  if val0_min < 0
    hl = line([-val0_min -val0_min], [0 max_h]);
    set(hl,'Color',[0 0.6667 1],'LineWidth',2);
    text(0.95*val0_max,0.85*max_h,'Max. observed value (inverse contrast) ',...
      'Color',[0 0.6667 1],'HorizontalAlignment','Right','FontSize',8)
  end
  
  % plot thresholds
  for j=1:n_alpha
    hl = line([val_th(n,j) val_th(n,j)], [0 max_h]);
    set(hl,'Color',col(j,:),'LineStyle','--');
    text(0.95*lim_x(2),(0.4+0.1*j)*max_h,['p<' num2str(alpha(j))],...
      'Color',col(j,:),'HorizontalAlignment','Right','FontSize',8)
  end
  
  ylabel('Frequency');
  xlabel(['Max ' name]);
  if corr
    title(['Distribution of maximum ' name],'FontWeight','bold');
  else
    title(['Distribution of ' name],'FontWeight','bold');
  end
  
  subplot(2,2,2*order)
  
  val_min = min(min(val_th(1:n,:)));
  val_max = max(max(val_th(1:n,:)));
  if val_max/val_min > 10
    hp = semilogy(1:n,val_th(1:n,:));
    yl = log10(ylim);
    ylim(10.^[floor(yl(1)) ceil(yl(2))])
  else
    hp = plot(1:n,val_th(1:n,:));
  end
  
  for j=1:n_alpha
    set(hp(j),'Color',col(j,:));
  end
  
  if corr
    title(['Corrected threshold of ' name],'FontWeight','bold')
  else
    title(['Uncorr. threshold of ' name],'FontWeight','bold')
  end
  ylabel('Threshold')
  xlabel('Permutations')   
end

%---------------------------------------------------------------
function [T, trRV] = calc_GLM(Y,xX,xCon,ind_mask,dim)
% compute T- or F-statistic using GLM
%
% Y        - masked data as vector
% xX       - design structure
% xCon     - contrast structure
% ind_mask - index of mask image
% dim      - image dimension
%
% Output:
% T        - T/F-values
% trRV     - df

c = xCon.c;

X = xX.W*xX.X;
pKX = pinv(X);

n_data = size(X,1);

Beta = Y*pKX';

res0 = Beta*(single(X'));
res0 = Y - res0; %-Residuals
res0 = res0.^2;
ResSS = double(sum(res0,2));
clear res0

trRV = n_data - rank(xX.X);
ResMS = ResSS/trRV;

%-Modify ResMS (a form of shrinkage) to avoid problems of very low variance
ResMS  = ResMS + 1e-3 * max(ResMS(isfinite(ResMS)));

T = zeros(dim);

if strcmp(xCon.STAT,'T')
  Bcov = pKX*pKX';
  con = Beta*c;

  T(ind_mask) = con./(eps+sqrt(ResMS*(c'*Bcov*c)));
else
  %-Compute ESS (ESS = Y'*Y - Beta'*X'*Y)
  h  = spm_FcUtil('Hsqr',xCon,xX.xKXs);
  ess = sum((h*Beta').^2,1)';
  MVM = ess/xCon.eidf;
  
  T(ind_mask) = MVM./ResMS;
end

%---------------------------------------------------------------
function [T, trRV] = calc_GLM_voxelwise(Y,xX,xC,xCon,ind_mask,dim,C,Pset,ind_X,pinv_method)
% compute T- or F-statistic using GLM
%
% Y        - masked data as vector
% xX       - design structure
% xC       - covariate structure
% xCon     - contrast structure
% ind_mask - index of mask image
% dim      - image dimension
%
% Output:
% T        - T/F-values
% trRV     - df

% use Matlab pinv by default
if nargin < 10, pinv_method = 1; end

c = xCon.c;
X = xX.W*xX.X;

n_data = size(X,1);
trRV = n_data - rank(xX.X);

T  = zeros(dim);
T0 = zeros(size(ind_mask));

% only permute C if this covariate is selected by the contrast
if ~isempty(Pset)
  found_ind_X = false;
  for j=1:numel(xC.cols)
    if ~isempty(find(ind_X==xC.cols(j), 1))
      found_ind_X = true;
    else
      found_ind_X = false;
    end
  end
  if found_ind_X
    C = C*full(Pset);
  end
end

% use pinv2 if it was faster
if pinv_method == 2
    pinvx  = @(x) pinv2(x); 
else
    pinvx  = @(x) pinv(x); 
end

% go through all voxels inside mask
for i=1:size(Y,1)

  % get X and replace defined columns with voxel-wise covariate
  Xi = X;
  for j=1:numel(xC.cols)
    ind = (X(:,xC.cols(j)) > 0);
    Xi(ind,xC.cols(j)) = C(i,ind)';
  end

  pKX = pinvx(Xi);

  Beta = pKX*Y(i,:)';

  res0 = Beta'*single(Xi)';
  res0 = res0 - Y(i,:); %-Residuals
  res0 = res0.^2;
  ResSS = double(sum(res0,2));
  
  ResMS = ResSS/trRV;
  %-Modify ResMS (a form of shrinkage) to avoid problems of very low variance
  ResMS  = ResMS + 1e-3 * max(ResMS(isfinite(ResMS)));
  
  if strcmp(xCon.STAT,'T')
    Bcov = pKX*pKX';
    con = Beta'*c;
  
    T0(i) = con./(eps+sqrt(ResMS*(c'*Bcov*c)));
  else
    error('F-test is not yet supported')
    %-Compute ESS
    h  = spm_FcUtil('Hsqr',xCon,xX.xKXs);
    ess = sum((h*Beta).^2,1)';
    MVM = ess/xCon.eidf;
    
    T(ind_mask) = MVM./ResMS;
  end
end

T(ind_mask) = T0;

%---------------------------------------------------------------
function xX = correct_xX(xX)
% sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design

% vector of covariates and nuisance variables
iCG = [xX.iC xX.iG];

% set columns with covariates and nuisance variables to zero
X = xX.X;
X(:,iCG) = 0;

ncol = size(X,2);

% calculate sum of columns
% The idea behind this is that for each factor the sum of all of its columns should be "1".
Xsum = zeros(size(X));
for i=1:ncol
  % only sum up columns without covariates and nuisance variables
  if isempty(find(iCG==i, 1))
    Xsum(:,i) = sum(X(:,1:i),2);
  end
end

% find columns where all entries are constant except zeros entries
% that indicate columns with covariates and nuisance variables
ind = find(any(diff(Xsum))==0 & sum(Xsum)>0);

% no more than 2 factors expected
if length(ind) > 2
  error('Weird design was found that cannot be analyzed correctly.');
end

% correction is only necessary if 2 factors (iH/iB) were found
if length(ind) > 1
  iF = cell(length(ind),1);

  j = 1;
  % skip columns with covariates and nuisance variables
  while find(iCG==j),  j = j + 1; end

  for i=j:length(ind)
    iF{i} = [j:ind(i)];
  
    j = ind(i)+1;
    % skip columns with covariates and nuisance variables
    while find(iCG==j), j = j + 1; end
  end
  
  % not sure whether this will always work but usually iB (subject effects) should be larger than iH (time effects)
%  if length(iF{1}) > length(iF{2})
if 0 % will be probably not always correct 
    xX.iB = iF{1};
    xX.iH = iF{2};
  else
    xX.iB = iF{2};
    xX.iH = iF{1};
  end
end

%---------------------------------------------------------------
function str = num2str_short(num)
% get shorther strings for continuous numbers with length > 4

if length(num) > 4
  % check whether vector consist of continuous numbers
  if all(diff(num)==1)
    str = [num2str(num(1)) ':' num2str(num(end))];
  else
    str = num2str(num);
  end
else
  str = num2str(num);
end

%---------------------------------------------------------------
function varargout = palm_moments(varargin)
% For a statistic that can be expressed as trace(A*W), for
% a sample size of n observations, this function returns the
% expected first three moments of the permutation distribution,
% without actually computing any permutation.
%
% [mu,sigsq,gamm1,gamm2] = palm_moments(G)
% [mu,sigsq,gamm1]       = palm_moments(A,W,n)
%
% Inputs:
% - G : A PxV array of observations of the G random variable.
%       The moments are unbiased and run along the 1st dimension.
%       Typical case is P = number of permutations and V = 
%       number of tests, e.g., voxels.
% - A : A square matrix for a multivariate proper statistic,
%       or a vector of A values for various univariate tests.
% - W : A square matrix for a multivariate proper statistic,
%       or a vector of W values for various univariate tests.
% - n : Sample size on which the statistic is based.
%
% Outputs:
% - mu    : Sample mean.
% - sigsq : Sample variance (unbiased).
% - gamm1 : Sample skewness (unbiased).
% - gamm2 : Sample kurtosis (unbiased).
%
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% For the estimators using trace(AW), the references are:
% * Kazi-Aoual F, Hitier S, Sabatier R, Lebreton J-D. Refined
%   approximations to permutation tests for multivariate
%   inference. Comput Stat Data Anal. 1995;20(94):643-656.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;4:497-511.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Mar/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin == 1
    
    % For a set of values of the random variable G, return the
    % first 4 moments.
    
    % Mean
    G     = varargin{1};
    n     = size(G,1);
    mu    = sum(G,1)/n;
    
    % Variance
    G0    = bsxfun(@minus,G,mu);
    ssq   = sum(G0.^2,1);
    sigsq = (ssq/(n-1));
    
    % Skewness
    s2    = ssq/n; % biased variance
    m3    = sum(G0.^3,1)/n;
    gamm1 = m3./s2.^1.5;
    gamm1 = gamm1 * sqrt((n-1)/n)*n/(n-2); % unbiased skewness

    % Kurtosis (normal dist = 3)
    if nargout == 4
        m4    = sum(G0.^4,1)/n;
        gamm2 = (m4./s2.^2);
        gamm2 = ((n+1)* gamm2 -3*(n-1))*(n-1)/((n-2)*(n-3))+3; % unbiased kurtosis
    else
        gamm2 = [];
    end
    
elseif nargin == 3
    
    % Compute the first three moments of the permutation distribution of
    % the statistic G = trace(AW), for n subjects, using the method in
    % Kazi-Aoual et al (1995). The variable names follow roughly the same
    % as in the paper.
    
    % Take inputs
    A = varargin{1};
    W = varargin{2};
    n = varargin{3};
    
    % If A and W are truly multivariate (i.e., square matrices), do as in
    % the original paper. Otherwise, make simplifications as these are all
    % scalars.
    if size(A,1) == size(A,2)
        
        % Some auxiliary variables for ET2:
        T    = trace(A);
        T_2  = trace(A^2);
        S_2  = sum(diag(A).^2);
        Ts   = trace(W);
        Ts_2 = trace(W^2);
        Ss_2 = sum(diag(W).^2);
        
        % Some auxiliary variables for ET3:
        T_3  = trace(A^3);
        S_3  = sum(diag(A).^3);
        U    = sum(A(:).^3);
        R    = diag(A)'*diag(A^2);
        B    = diag(A)'*A*diag(A);
        Ts_3 = trace(W^3);
        Ss_3 = sum(diag(W).^3);
        Us   = sum(W(:).^3);
        Rs   = diag(W)'*diag(W^2);
        Bs   = diag(W)'*W*diag(W);
        
    else
        
        % Some auxiliary variables for ET2:
        T    = A;
        T_2  = A.^2;
        S_2  = T_2;
        Ts   = W;
        Ts_2 = W.^2;
        Ss_2 = Ts_2;
        
        % Some auxiliary variables for ET3:
        T_3  = A.^3;
        S_3  = T_3;
        U    = T_3;
        R    = T_3;
        B    = T_3;
        Ts_3 = W.^3;
        Ss_3 = Ts_3;
        Us   = Ts_3;
        Rs   = Ts_3;
        Bs   = Ts_3;
    end
    
    % E(T):
    mu = T.*Ts/(n-1);
    
    % V(T):
    sigsq = 2*((n-1)*T_2-T.^2).*((n-1)*Ts_2-Ts.^2) / (n-1)^2/(n+1)/(n-2) ...
        + (n*(n+1)*S_2-(n-1)*(T.^2+2*T_2)) .* (n*(n+1)*Ss_2-(n-1)*(Ts.^2+2*Ts_2)) ...
        / (n+1)/n/(n-1)/(n-2)/(n-3);
    
    % E(T^3):
    ET3 = ...
        n^2*(n+1)*(n^2+15*n-4)*S_3.*Ss_3 ...
        + 4*(n^4-8*n^3+19*n^2-4*n-16)*U.*Us ...
        + 24*(n^2-n-4)*(U.*Bs+B.*Us) ...
        + 6*(n^4-8*n^3+21*n^2-6*n-24)*B.*Bs ...
        + 12*(n^4-n^3-8*n^2+36*n-48)*R.*Rs ...
        + 12*(n^3-2*n^2+9*n-12)*(T.*S_2.*Rs + R.*Ts.*Ss_2) ...
        + 3*(n^4-4*n^3-2*n^2+9*n-12)*T.*Ts.*S_2.*Ss_2 ...
        + 24*( (n^3-3*n^2-2*n+8)*(R.*Us+U.*Rs) ...
        + (n^3-2*n^2-3*n+12)*(R.*Bs+B.*Rs) ) ...
        + 12*(n^2-n+4)*(T.*S_2.*Us+U.*Ts.*Ss_2) ...
        + 6*(2*n^3-7*n^2-3*n+12)*(T.*S_2.*Bs+B.*Ts.*Ss_2) ...
        - 2*n*(n-1)*(n^2-n+4)*( (2*U+3*B).*Ss_3+(2*Us+3*Bs).*S_3 ) ...
        - 3*n*(n-1)^2*(n+4)*( (T.*S_2+4*R).*Ss_3+(Ts.*Ss_2+4*Rs).*S_3 ) ...
        + 2*n*(n-1)*(n-2)*( (T.^3+6*T.*T_2+8*T_3).*Ss_3 ...
        + (Ts.^3+6*Ts.*Ts_2+8*Ts_3).*S_3 ) ...
        + T.^3.*((n^3-9*n^2+23*n-14)*Ts.^3+6*(n-4).*Ts.*Ts_2+8*Ts_3) ...
        + 6*T.*T_2.*((n-4)*Ts.^3+(n^3-9*n^2+24*n-14)*Ts.*Ts_2+4*(n-3)*Ts_3) ...
        + 8*T_3.*(Ts.^3+3*(n-3).*Ts.*Ts_2+(n^3-9*n^2+26*n-22)*Ts_3) ...
        - 16*(T.^3.*Us+U.*Ts.^3)-6*(T.*T_2.*Us+U.*Ts.*Ts_2)*(2*n^2-10*n+16) ...
        - 8*(T_3.*Us+U.*Ts_3)*(3*n^2-15*n+16)-(T.^3.*Bs+B.*Ts.^3) ...
        * (6*n^2-30*n+24)-6*(T.*T_2.*Bs+B.*Ts.*Ts_2)*(4*n^2-20*n+24) ...
        - 8*(T_3.*Bs + B.*Ts_3)*(3*n^2-15*n+24) ...
        - (n-2)*( 24*(T.^3.*Rs+R.*Ts.^3)+6*(T.*T_2.*Rs+R.*Ts.*Ts_2)*(2*n^2-10*n+24) ...
        + 8*(T_3.*Rs+R.*Ts_3)*(3*n^2-15*n+24)+(3*n^2-15*n+6) ...
        .* (T.^3.*Ts.*Ss_2+T.*S_2.*Ts.^3) ...
        + 6*(T.*T_2.*Ts.*Ss_2+T.*S_2.*Ts.*Ts_2)*(n^2-5*n+6) ...
        + 48*(T_3.*Ts.*Ss_2+T.*S_2.*Ts_3) );
    ET3 = ET3/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5);
    
    % The coefficient "3" below is missing from Kazi-Aoual et al (1995), but it
    % is shown in the Supplementary Information of Minas and Montana (2014).
    gamm1 = (ET3 - 3*mu.*sigsq - mu.^3)./sigsq.^1.5;
    gamm2 = [];
else
    error('Incorrect number of arguments');
end

% Return results
varargout{1} = mu;
varargout{2} = sigsq;
varargout{3} = gamm1;
varargout{4} = gamm2;

%---------------------------------------------------------------
function pvals = palm_gamma(G,mu,sigsq,gamm1,rev,prepl)
% Return the p-values for a Gamma distribution, parameterised by
% its first three moments.
%
% pvals = palm_gamma(G,mu,s2,gamm1,rev)
% 
% Inputs:
% - G     : Statistics for which p-values are to be computed.
% - mu    : Distribution mean.
% - sigsq : Distribution standard deviation.
% - gamm1 : Distribution skewness.
% - rev   : Use if lower values of the statistic are evidence in
%           favour of the alternative.
% - prepl : Replacement for what otherwise would be zero p-values
%           in case of poor fits (e.g., statistic falls into the
%           part of the distribution that has pdf=0. In these cases
%           the p-value can be 1 or 1/(#perm) depending on which
%           tail and the sign of the skewness.
%
% Outputs:
% - pvals : p-values.
% 
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% Other references:
% * Mielke PW, Berry KJ, Brier GW. Application of Multi-Response
%   Permutation Procedures for Examining Seasonal Changes in
%   Monthly Mean Sea-Level Pressure Patterns. Mon Weather Rev.
%   1981;109(1):120-126.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;7(6):450-470.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Note that there are no argument checking for speed, but
% sizes of all inputs need to be the same, or the moments need to
% be all scalars.

if gamm1 == 0
    
    % If not skewed, use a normal approximation.
    G     = (G - mu)./sigsq.^.5;
    pvals = erfc(G/sqrt(2))/2;
    
else
    
    % Standardise G, so that all becomes a function of the skewness.
    G    = (G - mu)./sigsq.^.5;
    
    % Gamma distribution parameters (Minas & Montana, 2014).
    kpar = 4/gamm1.^2;
    tpar = gamm1/2;
    cpar = -2/gamm1;
     
    % Actual p-value. If there are negatives here, the probability can
    % have an imaginary part, which is dealt with later.
    if rev
        if gamm1 > 0
            tail = 'lower';
        else
            tail = 'upper';
        end
    else
        if gamm1 > 0
            tail = 'upper';
        else
            tail = 'lower';
        end
    end
    pvals = gammainc((G-cpar)./tpar,kpar,tail);
    
    % Deal with imaginary parts.
    if ~ isreal(pvals)
        iidx = imag(pvals) ~= 0;
        if rev
            if gamm1 > 0
                pvals(iidx) = prepl;
            else
                pvals(iidx) = 1;
            end
        else
            if gamm1 > 0
                pvals(iidx) = 1;
            else
                pvals(iidx) = prepl;
            end
        end
    end
end

%---------------------------------------------------------------
function Z = palm_gtoz(G,df1,df2)
% Convert a G-statistic (or any of its particular cases)
% to a z-statistic (normally distributed).
%
% Usage:
% Z = palm_gtoz(G,df1,df2)
%
% Inputs:
% G        : G statistic.
% df1, df2 : Degrees of freedom (non-infinite).
% 
% Outputs:
% Z        : Z-score
%
% If df2 = NaN and df1 = 1, G is treated as Pearson's r.
% If df2 = NaN and df1 > 1, G is treated as R^2.
% If df2 = NaN and df1 = 0, G is treated as z already.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Jan/2014
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Note that for speed, there's no argument checking.

% If df2 is NaN, this is r, R^2, or z already
if isnan(df2(1))
    
  if df1 == 0
      
    % If df1 is zero, this is already a z-stat (this is here more for
    % compatibility).
    Z = G;
      
  elseif df1 == 1
      
    % If rank(C) = 1, i.e., df1 = 1, this is r, so
    % do a Fisher's r-to-z stransformation
    Z = atanh(G);
      
  elseif df1 > 1
      
    % If rank(C) > 1, i.e., df1 > 1, this is R^2, so
    % use a probit transformation.
    Z = -erfcinv(2*G)*sqrt(2); %Z = norminv(G);
      
  end

else
  siz = size(G);
  Z = zeros(siz);
  df2 = bsxfun(@times,ones(siz),df2);
  if df1 == 1
    
    % Deal with precision issues working on each
    % tail separately
    idx = G > 0;
    Z( idx) =  erfcinv(2*palm_gcdf(-G( idx),1,df2( idx)))*sqrt(2);
    Z(~idx) = -erfcinv(2*palm_gcdf( G(~idx),1,df2(~idx)))*sqrt(2);
    
  elseif df1 == 0
    
    % If df1 is zero, this is already a z-stat (this is here more for
    % compatibility).
    Z = G;
    
  else
    
    % G-vals above the upper half are treated as
    % "upper tail"; otherwise, "lower tail".
    thr = (1./betainv(.5,df2/2,df1/2)-1).*df2/df1;
    idx = G > thr;
    
    % Convert G-distributed variables to Beta-distributed
    % variables with parameters a=df1/2 and b=df2/2
    B = (df1.*G./df2)./(1+df1.*G./df2);
    a = df1/2;
    b = df2/2;
    
    % Convert to Z through a Beta incomplete function
    %Z( idx) = -erfinv(2*betainc(1-B( idx),b( idx),a)-1)*sqrt(2);
    %Z(~idx) =  erfinv(2*betainc(  B(~idx),a,b(~idx))-1)*sqrt(2);
    Z( idx) =  erfcinv(2*betainc(1-B( idx),b( idx),a))*sqrt(2);
    Z(~idx) = -erfcinv(2*betainc(  B(~idx),a,b(~idx)))*sqrt(2);
    
  end
end

%---------------------------------------------------------------
% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2012 Kurt Hornik
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {} betainv (@var{x}, @var{a}, @var{b})
% For each element of @var{x}, compute the quantile (the inverse of
% the CDF) at @var{x} of the Beta distribution with parameters @var{a}
% and @var{b}.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: Quantile function of the Beta distribution

function inv = betainv (x, a, b)

if (nargin ~= 3)
  print_usage ();
end

if (~isscalar (a) || ~isscalar (b))
  [retval, x, a, b] = common_size (x, a, b);
  if (retval > 0)
    error ('betainv: X, A, and B must be of common size or scalars');
  end
end

if (iscomplex (x) || iscomplex (a) || iscomplex (b))
  error ('betainv: X, A, and B must not be complex');
end

if (isa (x, 'single') || isa (a, 'single') || isa (b, 'single'))
  inv = zeros (size (x), 'single');
else
  inv = zeros (size (x));
end

k = (x < 0) | (x > 1) | ~(a > 0) | ~(b > 0) | isnan (x);
inv(k) = NaN;

k = (x == 1) & (a > 0) & (b > 0);
inv(k) = 1;

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
  if (~isscalar (a) || ~isscalar (b))
    a = a(k);
    b = b(k);
    y = a ./ (a + b);
  else
    y = a / (a + b) * ones (size (k));
  end
  x = x(k);
  
  if (isa (y, 'single'))
    myeps = eps ('single');
  else
    myeps = eps;
  end
  
  l = find (y < myeps);
  if (any (l))
    y(l) = sqrt (myeps) * ones (length (l), 1);
  end
  l = find (y > 1 - myeps);
  if (any (l))
    y(l) = 1 - sqrt (myeps) * ones (length (l), 1);
  end
  
  y_old = y;
  for i = 1 : 10000
    h     = (betacdf (y_old, a, b) - x) ./ betapdf (y_old, a, b);
    y_new = y_old - h;
    ind   = find (y_new <= myeps);
    if (any (ind))
      y_new (ind) = y_old (ind) / 10;
    end
    ind = find (y_new >= 1 - myeps);
    if (any (ind))
      y_new (ind) = 1 - (1 - y_old (ind)) / 10;
    end
    h = y_old - y_new;
    if (max (abs (h)) < sqrt (myeps))
      break;
    end
    y_old = y_new;
  end
  
  inv(k) = y_new;
end


%---------------------------------------------------------------
function gcdf = palm_gcdf(G,df1,df2)
% Convert a pivotal statistic computed with 'pivotal.m'
% (or simplifications) to a parametric p-value.
% The output is 1-p, i.e. the CDF.
% 
% Usage:
% cdf = palm_gcdf(G,df1,df2)
% 
% Inputs:
% G        : G or Z statistic.
% df1, df2 : Degrees of freedom (non infinite).
%            df1 must be a scalar
%            For z, use df1 = 0.
%            For Chi2, use df1 = -1, and df2 as the df.
% 
% Outputs:
% cdf      : Parametric cdf (1-p), based on a
%            t, F, z or Chi2 distribution.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Note that for speed, there's no argument checking,
% and some lines are repeated inside the conditions.

if df1 > 1
  
  % G or F
  df2 = bsxfun(@times,ones(size(G)),df2);
  B = (df1.*G./df2)./(1+df1.*G./df2);
  gcdf = betainc(B,df1/2,df2/2);

elseif df1 == 1
  
  % Student's t, Aspin's v
  df2 = bsxfun(@times,ones(size(G)),df2);
  ic = df2 == 1;
  in = df2 > 1e7;
  ig = ~(ic|in);
  gcdf = zeros(size(G));
  if any(ig(:))
    gcdf(ig) = betainc(1./(1+G(ig).^2./df2(ig)),df2(ig)/2,.5)/2;
  end
  ig = G > 0 & ig;
  gcdf(ig) = 1 - gcdf(ig);
  if any(ic(:))
    gcdf(ic) = .5 + atan(G(ic))/pi;
  end
  if any(in(:))
    gcdf(ic) = palm_gcdf(G(in),0);
  end

elseif df1 == 0
  
  % Normal distribution
  gcdf = erfc(-G/sqrt(2))/2;
  
elseif df1 < 0
  
  % Chi^2, via lower Gamma incomplete for precision and speed
  %df2 = bsxfun(@times,ones(size(G)),df2);
  gcdf = gammainc(G/2,df2/2,'lower');
  
end

%---------------------------------------------------------------
% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2012 Kurt Hornik
% Copyright (C) 2010 Christos Dimitrakakis
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {} betapdf (@var{x}, @var{a}, @var{b})
% For each element of @var{x}, compute the probability density function (PDF)
% at @var{x} of the Beta distribution with parameters @var{a} and @var{b}.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>, CD <christos.dimitrakakis@gmail.com>
% Description: PDF of the Beta distribution

function pdf = betapdf (x, a, b)

if (nargin ~= 3)
  print_usage ();
end

if (~isscalar (a) || ~isscalar (b))
  [retval, x, a, b] = common_size (x, a, b);
  if (retval > 0)
    error ('betapdf: X, A, and B must be of common size or scalars');
  end
end

if (iscomplex (x) || iscomplex (a) || iscomplex (b))
  error ('betapdf: X, A, and B must not be complex');
end

if (isa (x, 'single') || isa (a, 'single') || isa (b, 'single'))
  pdf = zeros (size (x), 'single');
else
  pdf = zeros (size (x));
end

k = ~(a > 0) | ~(b > 0) | isnan (x);
pdf(k) = NaN;

k = (x > 0) & (x < 1) & (a > 0) & (b > 0) & ((a ~= 1) | (b ~= 1));
if (isscalar (a) && isscalar (b))
  pdf(k) = exp ((a - 1) * log (x(k))...
                + (b - 1) * log (1 - x(k))...
                + lgamma (a + b) - lgamma (a) - lgamma (b));
else
  pdf(k) = exp ((a(k) - 1) .* log (x(k))...
                + (b(k) - 1) .* log (1 - x(k))...
                + lgamma (a(k) + b(k)) - lgamma (a(k)) - lgamma (b(k)));
end

% Most important special cases when the density is finite.
k = (x == 0) & (a == 1) & (b > 0) & (b ~= 1);
if (isscalar (a) && isscalar (b))
  pdf(k) = exp (lgamma (a + b) - lgamma (a) - lgamma (b));
else
  pdf(k) = exp (lgamma (a(k) + b(k)) - lgamma (a(k)) - lgamma (b(k)));
end

k = (x == 1) & (b == 1) & (a > 0) & (a ~= 1);
if (isscalar (a) && isscalar (b))
  pdf(k) = exp (lgamma (a + b) - lgamma (a) - lgamma (b));
else
  pdf(k) = exp (lgamma (a(k) + b(k)) - lgamma (a(k)) - lgamma (b(k)));
end

k = (x >= 0) & (x <= 1) & (a == 1) & (b == 1);
pdf(k) = 1;

% Other special case when the density at the boundary is infinite.
k = (x == 0) & (a < 1);
pdf(k) = Inf;

k = (x == 1) & (b < 1);
pdf(k) = Inf;

%---------------------------------------------------------------
function [errorcode, varargout] = common_size (varargin)

if (nargin < 2)
  error ('common_size: only makes sense if nargin >= 2');
end

% Find scalar args.
nscal = cellfun (@numel, varargin) ~= 1;

i = find (nscal, 1);

if (isempty (i))
  errorcode = 0;
  varargout = varargin;
else
  match = cellfun (@size_equal, varargin, repmat(varargin(i),size(varargin)));
  if (any (nscal & ~match))
    errorcode = 1;
    varargout = varargin;
  else
    errorcode = 0;
    if (nargout > 1)
      scal = ~nscal;
      varargout = varargin;
      dims = size (varargin{i});
      for s = find(scal)
        varargout{s} = repmat(varargin{s}, dims);
      end
    end
  end
end

%---------------------------------------------------------------
% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2012 Kurt Hornik
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {} betacdf (@var{x}, @var{a}, @var{b})
% For each element of @var{x}, compute the cumulative distribution function
% (CDF) at @var{x} of the Beta distribution with parameters @var{a} and
% @var{b}.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: CDF of the Beta distribution

function cdf = betacdf (x, a, b)

if (nargin ~= 3)
  print_usage ();
end

if (~isscalar (a) || ~isscalar (b))
  [retval, x, a, b] = common_size (x, a, b);
  if (retval > 0)
    error ('betacdf: X, A, and B must be of common size or scalars');
  end
end

if (iscomplex (x) || iscomplex (a) || iscomplex (b))
  error ('betacdf: X, A, and B must not be complex');
end

if (isa (x, 'single') || isa (a, 'single') || isa (b, 'single'))
  cdf = zeros (size (x), 'single');
else
  cdf = zeros (size (x));
end

k = isnan (x) | ~(a > 0) | ~(b > 0);
cdf(k) = NaN;

k = (x >= 1) & (a > 0) & (b > 0);
cdf(k) = 1;

k = (x > 0) & (x < 1) & (a > 0) & (b > 0);
if (isscalar (a) && isscalar (b))
  cdf(k) = betainc (x(k), a, b);
else
  cdf(k) = betainc (x(k), a(k), b(k));
end

%---------------------------------------------------------------
function eq = size_equal(a,b)
eq = isequal(size(a),size(b));

%---------------------------------------------------------------
function a = iscomplex(X)
a = ~isreal(X);

%---------------------------------------------------------------
function varargout = lgamma(varargin)
varargout{1:nargout} = gammaln(varargin{:});

%---------------------------------------------------------------
function X = pinv2(A,tol)
%PINV2   Pseudoinverse.
%   X = PINV2(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%
%   PINV2(A,TOL) treats all singular values of A that are less than TOL as
%   zero. By default, TOL = max(size(A)) * eps(norm(A)).
%
% This function is a modified version of pinv from Matlab but uses a faster
% svdecon function from Vipin Vijayan.

[U,S,V] = svdecon(A);
s = diag(S);
if nargin < 2 
    tol = max(size(A)) * eps(norm(s,inf));
end
r1 = sum(s > tol)+1;
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
s = 1./s(:);
X = (V.*s.')*U';


function [U,S,V] = svdecon(X)
% Input:
% X : m x n matrix
%
% Output:
% X = U*S*V'
%
% Description:
% Does equivalent to svd(X,'econ') but faster
%
% Vipin Vijayan (2014)

[m,n] = size(X);
if  m <= n
    C = X*X';
    [U,D] = eig(C);
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    U = U(:,ix);    
    
    if nargout > 2
        V = X'*U;
        s = sqrt(d);
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
    C = X'*X; 
    [V,D] = eig(C);
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    V = V(:,ix);    
    
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(d);
    U = bsxfun(@(x,c)x./c, U, s');
    S = diag(s);
end
