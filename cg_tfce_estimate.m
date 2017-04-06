function cg_tfce_estimate(job)
% main TFCE function for estimating TFCE statistics
%
% FORMAT cg_tfce_estimate(job)
% job - job from tbx_cfg_tfce_estimate
% 
%_______________________________________________________________________
% Christian Gaser
% $Id$

% convert to z-statistic
convert_to_z = 0;

% display permuted design matrix (otherwise show t distribution)
show_permuted_designmatrix = 1;

% allow to test permutations without analyzing data
test_mode = 0;

% define stepsize for tfce
n_steps_tfce = 100;
    
% colors and alpha levels
col   = [0.25 0 0; 1 0 0; 1 0.75 0];
alpha = [0.05 0.01 0.001];
    
% give same results each time
if exist('rng','builtin') == 5
  rng('default')
  rng(0)
else
  rand('state',0);
end

% tolerance for comparing real numbers
tol = 1e-4;
    
% check spm version
if strcmp(spm('ver'),'SPM12')
  spm12 = 1;
else
  spm12 = 0;
end

load(job.spmmat{1});
cwd = fileparts(job.spmmat{1});

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
try
  xCon0 = SPM.xCon(Ic0(1));
catch
  [Ic0,xCon] = spm_conman(SPM,'T&F',Inf,...
        '  Select contrast(s)...',' ',1);
  SPM.xCon = xCon;
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

% find exchangeability block labels for longitudinal designs (paired t-test, flexible factorial)
repeated_anova = ~isempty(xX.iB);
if repeated_anova
  [rw,cl] = find(xX.I == length(xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_subj)
  exch_block_labels = xX.I(:,cl(1));     % column from above contains the subject factor

  % check that labels are defined for each block
  for i=1:n_data
    groupListed(exch_block_labels(i)) = true;
  end
  for i = 1: max(exch_block_labels)
    if ~groupListed(i) 
      fprintf('Error: block %d must be assigned to at least one design row in the blocks file.\n',i);
      return
    end
  end
  fprintf('Please note that permutation is only done within subjects for repeated Anova.\n',i);
else
  exch_block_labels = ones(1,n_data);
end

% check for meshes
if spm12
  if spm_mesh_detect(VY)
    mesh_detected = 1;
  else
    mesh_detected = 0;
  end
else
  mesh_detected = 0;
end

% set E according to type data
if job.tbss || mesh_detected
  E = 1.0; H = 2.0;
else
  E = 0.5; H = 2.0;
end

if ~test_mode
  % check for mask image that should exist for any analysis
  if exist(fullfile(cwd, 'mask.img'))
    file_ext = '.img';
  elseif exist(fullfile(cwd, 'mask.nii'))
    file_ext = '.nii';
  elseif exist(fullfile(cwd, 'mask.gii'))
    file_ext = '.gii';
  else
    spm('alert!',sprintf('WARNING: No mask file found. Switch to test mode.\n\n'),'',spm('CmdLine'),1);
    test_mode = 1;
  end
end
  
if ~test_mode
  % get mask file
  if isempty(job.mask)
    maskname = fullfile(cwd,['mask' file_ext]);
  else
    maskname = job.mask{1};
  end
      
  % load mask
  try
    if spm12
      Vmask = spm_data_hdr_read(maskname);
    else
      Vmask = spm_vol(maskname);
    end
  catch
    if mesh_detected
      maskname = spm_select(1,'mesh','select mask image');
    else
      maskname = spm_select(1,'image','select mask image');
    end
    if spm12
      Vmask = spm_data_hdr_read(maskname);
    else
      Vmask = spm_vol(maskname);
    end
  end
    
  % if first image was not found you have to select all files again
  if ~exist(VY(1).fname);
  
    n = size(SPM.xY.VY,1);
    if mesh_detected
      P = spm_select(n,'mesh','select images');
    else
      P = spm_select(n,'image','select images');
    end
    
    if spm12
      VY = spm_data_hdr_read(P);
    else
      VY = spm_vol(P);
    end
    
    %-Apply gSF to memory-mapped scalefactors to implement scaling
    %--------------------------------------------------------------------------
    for i = 1:n
      VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*SPM.xGX.gSF(i); % FIXME % for meshes
    end
    
    SPM.xY.VY = VY;
      
    % update SPM
    if size(SPM.xY.VY,1)==n
      save(job.spmmat{1},'SPM');
    else
      fprintf('ERROR: Number of files is not correct\n');
      return
    end
  end
        
  % check whether mask images fits to the data
  if mesh_detected, dim_index = 1; else dim_index=1:3; end
  if sum(sum((Vmask.mat-VY(1).mat).^2)) > 1e-6 || any(Vmask.dim(dim_index) ~= VY(1).dim(dim_index))
    error('Mask must have the same dimensions and orientation as the data.');
  end

  % read mask and data
  if spm12
    mask = spm_data_read(Vmask);
  else
    mask = spm_read_vols(Vmask);
  end
  
  ind_mask = find(mask>0);
  n = numel(VY);

  if ~isempty(ind_mask)
    Y = zeros([length(ind_mask) n],'single');
  
    % load data
    for i=1:n
      if spm12
        tmp = spm_data_read(VY(i));
      else 
        tmp = spm_read_vols(VY(i));
      end
      Y(:,i) = tmp(ind_mask);
    end
    clear tmp;
  
    % whitening matrix
    W = single(full(xX.W));
    Y = Y*W;
  else
    error('Empty mask.');
  end

  t0 = zeros(Vmask.dim);
  t  = zeros(Vmask.dim);

end % if ~test_mode

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
  
  % get contrast and name
  c = xCon.c;
  
  [indi, indj] = find(c~=0);
  n_contrasts_lines = size(c,2);
  ind_X = unique(indi)';
  
  % check for contrasts that are defined for columns with subject effects
  if ~isempty(xX.iB)
    if max(ind_X) > min(xX.iB)
      fprintf('ERROR: No contrasts on subjects/block effects allowed.\n');
      return
    end
  end

  c_name  = deblank(xCon.name);

  % find exchangeability blocks using contrasts without zero values
  exch_blocks   = c(ind_X);
  
  n_exch_blocks = length(ind_X);
  
  % check for exchangeability blocks and design matrix
  if n_exch_blocks == 1
    n_cond = length(find(xX.iH==ind_X)); % check whether the contrast is defined at columns for condition effects
  else
    n_cond = 0;
    n_data_cond = [];
    for k=1:length(xX.iH)
      n_data_cond = [n_data_cond sum(xX.X(:,xX.iH(k)))];
    end
    for j=1:n_exch_blocks
      c_exch_blocks = find(c==exch_blocks(j));
      for k=1:length(c_exch_blocks)
        n_cond = n_cond + length(find(xX.iH==c_exch_blocks(k)));
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
      elseif sum(n_data_cond(find(c==exch_blocks(1)))) == sum(n_data_cond(find(c==exch_blocks(2))))
        use_half_permutations = 1;
        fprintf('Equal sample sizes: half of permutations are used.\n');
      end
    end
  end
  
  ind_exch_blocks = cell(n_exch_blocks,1);
  for j=1:n_exch_blocks
    if strcmp(xCon.STAT,'T')
      ind_exch_blocks{j} = find(c==exch_blocks(j));
    else
      ind_exch_blocks{j} = ind_X(j);
    end
  end
      
  % check design
  switch n_cond
  case 0 % correlation
    if n_exch_blocks == 2
      disp('Interaction design between two regressors found. This should work, but is not yet fully tested and I am unsure whether the # of max. permutations is correctly estimated.')
    else
      if repeated_anova
        fprintf('Repeated Anova with contrast for covariate found\n');
      else
        fprintf('Multiple regression design found\n');
      end
    end
    label = 1:n_data;
  case 1 % one-sample t-test
    fprintf('One sample t-test found\n');
    % use exchangeability blocks for labels
    label = zeros(1,n_data);
    for j=1:n_exch_blocks
      for k=1:length(ind_exch_blocks{j})
        label(find(xX.X(:,ind_exch_blocks{j}(k))==1)) = j;
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
        label(find(xX.X(:,ind_exch_blocks{j}(k))==1)) = j;
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
      if (n_cond == 2) & (single_subject == 1)
        n_perm_full = n_data_with_contrast;
      else
        n_perm_full = realmax;
      end
    end
    
    % find where data are defined for that contrast
    if ~isempty(find(xX.iH == ind_X(1)))
      % first checking whether contrasts are defined for iH
      ind_data_defined = find(any(xX.X(:,xX.iH(ind_X)),2));
    else
      % if not check iC
      if length(ind_X) == 1 % for one contrast covariate column use all subjects
        ind_data_defined = 1:size(xX.X,1);
      else % check wehere contrast is defined
        ind_data_defined = find(any(xX.X(:,ind_X),2));
      end
    end
    
    % and restrict exchangeability block labels to those rows
    exch_block_labels_new = exch_block_labels(ind_data_defined);

    % Repated Anova: n_perm = n_cond1!*n_cond2!*...*n_condk!
    % for a full model where each condition is defined for all subjects the easier
    % estimation is: n_perm = (n_cond!)^n_subj
    % check that no regression analysis inside repeated anova is used
    if repeated_anova & n_cond~=0
      n_subj = max(exch_block_labels_new);
      n_perm_full = 1;
      for k=1:n_subj
        n_cond_subj = length(find(exch_block_labels_new == k));
        n_perm_full = n_perm_full*factorial(n_cond_subj);
      end
    else
      n_perm_full = round(n_perm_full);
    end
    
  else  % one-sample t-test: n_perm = 2^n
    n_perm_full = 2^n_data_with_contrast;
    exch_block_labels_new = exch_block_labels;
  end

  fprintf('# full permutations: %d\n',n_perm_full);
  if use_half_permutations
    fprintf('Equal sample size found: Use half permutations.\n');
  end
  fprintf('Exchangeability blocks: ');
  fprintf('%d ',unique(cell2mat(ind_exch_blocks)));
  fprintf('\n');
  fprintf('# of conditions: %d\n',n_cond);
  
  n_perm = min([n_perm n_perm_full]);
         
  % deal with interaction design
  if n_cond == 0 & n_exch_blocks == 2
    ind_X2 = [ind_X xX.iH];
  else
    ind_X2 = ind_X;
  end
  
  % Guttman partioning of design matrix into effects of interest X and nuisance variables Z
  X = xX.X(:,ind_X2);
  ind_Z = [xX.iH xX.iC xX.iB xX.iG];
  ind_Z(ind_X2) = [];
  Z = xX.X(:,ind_Z);
  
  Hz = Z*pinv(Z);
  Rz = eye(size(X,1)) - Hz;
  
  % if Hz is zero or Ic is empty then no confounds were found and we can skip the time-consuming
  % Freedman-Lane permutation
  if (all(~any(Hz)) | isempty(xX.iC)) & job.freedman_lane
    fprintf('No nuisance variables were found: Use Draper-Stoneman permutation.\n\n');
    job.freedman_lane = 0;
  end
  
  if ~test_mode

    if job.freedman_lane
      str_permutation_method = 'Freedman-Lane';
    else
      str_permutation_method = 'Draper-Stoneman';
    end

    % compute unpermuted t/F-map
    [t0, df2] = calc_GLM(Y,xX,xCon,ind_mask,VY(1).dim);
    
    mask_0 = (t0 == 0);
    mask_1 = (t0 ~= 0);
    mask_P = (t0 > 0);
    mask_N = (t0 < 0);

    df1 = size(xCon.c,2);

    %---------------------------------------------------------------
    % save unpermuted map
    %---------------------------------------------------------------

    % prepare output files
    Vt = VY(1);
    Vt.dt(1)    = 16;
    Vt.pinfo(1) = 1;

    name = sprintf('%s_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d %s',xCon.STAT,Ic, str_permutation_method);
    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,t0);
    else
      spm_write_vol(Vt,t0);
    end

    % transform to z statistic
    if convert_to_z
      % use faster z-transformation of SPM for T-statistics
      if strcmp(xCon.STAT,'T')
        t0 = spm_t2z(t0,df2);
      else
        t0 = palm_gtoz(t0,df1,df2);
      end
    end
    
    % remove all NaN and Inf's
    t0(isinf(t0) | isnan(t0)) = 0;

    % sometimes z-transformation produces neg. values even for F-statistics
    if strcmp(xCon.STAT,'F')
      t0(find(t0 < 0)) = 0;
    end

    % get dh for unpermuted map
    dh = max(abs(t0(:)))/n_steps_tfce;
  
    % calculate tfce of unpermuted t-map
    if mesh_detected
      if ~isstruct(SPM.xVol.G)
        SPM.xVol.G = gifti(SPM.xVol.G);
      end
      tfce0 = tfce_mesh(SPM.xVol.G.faces, t0, dh)*dh;
    else
      % only estimate neg. tfce values for non-positive t-values
      if min(t0(:)) < 0
        tfce0 = tfceMex_pthread(t0,dh,E,H,1,job.singlethreaded)*dh;
      else
        tfce0 = tfceMex_pthread(t0,dh,E,H,0,job.singlethreaded)*dh;
      end
    end
  
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
    t_th         = [];
    tfce_min     = [];
    tfce_max     = [];
    tfce_max_th  = [];
    tfce_th      = [];
  
  end % test_mode
  
  label_matrix = [];
  
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
  
    text(0.5,0.25,spm_str_manip(SPM.swd,'a50'),...
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
  if ~test_mode, cg_progress('Init',n_perm,'Calculating','Permutations'); end
  
  % update interval for progress bar
  progress_step = max([1 round(n_perm/100)]);

  ind_label_gt0 = find(label > 0);
  unique_labels = unique(label);
  n_unique_labels = length(unique_labels);
  
  i = 1;
  while i<=n_perm
  
    % randomize subject vector
    if i==1 % first permutation is always unpermuted model
      if n_cond == 1 % one-sample t-test
        rand_label = ones(1,n_data_with_contrast);
      else % correlation or Anova
        rand_order = ind_label;
        rand_order_sorted = rand_order;
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
        for k = 1:max(exch_block_labels_new)
          ind_block   = find(exch_block_labels_new == k);
          n_per_block = length(ind_block);
          rand_order(ind_block) = ind_label(ind_block(randperm(n_per_block)));
        end
        
        % go through defined labels and sort inside
        for k=1:n_unique_labels
          ind_block = find(label(ind_label_gt0) == unique_labels(k));
          rand_order_sorted(ind_block) = sort(rand_order(ind_block));
        end
        
        % check whether this permutation was already used
        while any(ismember(label_matrix,rand_order_sorted,'rows'))
          for k = 1:max(exch_block_labels_new)
            ind_block   = find(exch_block_labels_new == k);
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
        Pset(k,k) = rand_label(k);
      end
    else % correlation or Anova
      for k=1:n_data_with_contrast
        Pset(rand_order_sorted(k),ind_label(k)) = 1;
      end
    end

    % add Stop button after 20 iterations
    try % use try commands to allow batch mode without graphical output
      if i==21
        hStopButton = uicontrol(Fgraph,...
          'position',[10 10 70 20],...
          'style','toggle',...
          'string','Stop',...
          'backgroundcolor',[1 .5 .5]); % light-red
      end
    
      if i>=21
        stopStatus = get(hStopButton,'value');
      end
    
      % check Stop status
    if (stopStatus == true)
        fprintf('Stopped after %d iterations.\n',i);
        break; % stop the permutation loop
      end
    end
      
    % change design matrix according to permutation order
    % only permute columns, where contrast is defined
    Xperm = xX.X;
    
    % Draper-Stone permutation of design matrix only
    if ~job.freedman_lane
      Xperm(:,ind_X2) = Pset*Xperm(:,ind_X2);
    end

    Xperm_debug = xX.X;
    Xperm_debug(:,ind_X2) = Pset*Xperm_debug(:,ind_X2);

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
      
      % use different colors for indicated columns
      Xperm_debug(:,xX.iH) = 16*Xperm_debug(:,xX.iH);
      Xperm_debug(:,xX.iC) = 24*Xperm_debug(:,xX.iC);
      Xperm_debug(:,xX.iB) = 32*Xperm_debug(:,xX.iB);
      Xperm_debug(:,xX.iG) = 48*Xperm_debug(:,xX.iG);

      if n_cond==1 % one-sample t-test
        for j=1:n_data_with_contrast
          if rand_label(j) > 0
            Xperm_debug(ind_label(j),ind_X2) = 60*rand_label(j)*Xperm_debug(ind_label(j),ind_X2);
          else
            Xperm_debug(ind_label(j),ind_X2) = 56*rand_label(j)*Xperm_debug(ind_label(j),ind_X2);
          end
        end
      else % correlation or Anova
        % scale exchangeability blocks also to values 0.8..1
        val = Xperm_debug(:,ind_X2);
        ind0 = find(val==0);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
      
        % rescue zero entries
        val(ind0) = 0;
      
        Xperm_debug(:,ind_X2) = 60*val;
      end

    end
          
    
    % display permuted design matrix
    try
      if show_permuted_designmatrix
        subplot(2,2,3);
        image(Xperm_debug); axis off
        title('Permuted design matrix','FontWeight','bold');
      
        % use different colormap for permuted design matrix
        cmap = jet(64);
      
        % zero values should be always black
        cmap(1,:) = [0 0 0];
        colormap(cmap)
      
        subplot(2,2,4); axis off
      
        % color-coded legend
        y = 1.0;
        text(-0.2,y, 'Columns of design matrix: ', 'Color',cmap(1, :),'FontWeight','Bold','FontSize',10); y = y - 0.10;
        text(-0.2,y,['Exchangeability blocks: ' num2str_short(unique(cell2mat(ind_exch_blocks))')], 'Color',cmap(60,:),'FontWeight','Bold','FontSize',10); y = y - 0.05;
        if ~isempty(xX.iH)
          text(-0.2,y, ['iH - Indicator variables: ' num2str_short(xX.iH)], 'Color',cmap(16,:),'FontWeight','Bold','FontSize',10);
          y = y - 0.05; 
        end
        if ~isempty(xX.iC)
          text(-0.2,y, ['iC - Covariates: ' num2str_short(xX.iC)], 'Color',cmap(24,:),'FontWeight','Bold','FontSize',10);
          y = y - 0.05;
        end
        if ~isempty(xX.iB)
          text(-0.2,y, ['iB - Block effects: ' num2str_short(xX.iB)], 'Color',cmap(32,:),'FontWeight','Bold','FontSize',10);
          y = y - 0.05;
        end
        if ~isempty(xX.iG)
          text(-0.2,y, ['iG - Nuisance variables: ' num2str_short(xX.iG)], 'Color',cmap(48,:),'FontWeight','Bold','FontSize',10);
          y = y - 0.05;
        end
      end
    end
    
    if ~test_mode
      % calculate permuted t-map
      if i==1
        t    = t0;
        tfce = tfce0;
      else
        xXperm   = xX;
        xXperm.X = Xperm;

        % Freedman-Lane permutation of data
        if job.freedman_lane
          t = calc_GLM(Y*(Pset'*Rz+Hz),xXperm,xCon,ind_mask,VY(1).dim);
        else
          t = calc_GLM(Y,xXperm,xCon,ind_mask,VY(1).dim);
        end
        
        if convert_to_z
          % use faster z-transformation of SPM for T-statistics
          if strcmp(xCon.STAT,'T')
            t(mask_1) = spm_t2z(t(mask_1),df2);
          else
            t(mask_1) = palm_gtoz(t(mask_1),df1,df2);
          end
        end
      
        % remove all NaN and Inf's
        t(isinf(t) | isnan(t)) = 0;

        % use individual dh
        dh = max(abs(t(:)))/n_steps_tfce;
        
        % compute tfce
        if mesh_detected
          tfce = tfce_mesh(SPM.xVol.G.faces, t, dh)*dh;
        else
          % only estimate neg. tfce values for non-positive t-values
          if min(t(:)) < 0
            tfce = tfceMex_pthread(t,dh,E,H,1,job.singlethreaded)*dh;
          else
            tfce = tfceMex_pthread(t,dh,E,H,0,job.singlethreaded)*dh;
          end
        end
        
      end
    end % test_mode
    
    % update label_matrix for checking of unique permutations
    if use_half_permutations
      label_matrix = [label_matrix; rand_order_sorted; [rand_order_sorted(find(label(ind_label) == 2)) rand_order_sorted(find(label(ind_label) == 1))]];

      if ~test_mode
        % maximum statistic
        tfce_max = [tfce_max max(tfce(:)) -min(tfce(:))];
        t_max    = [t_max max(t(:)) -min(t(:))];
        tfce_min = [tfce_min min(tfce(:)) -max(tfce(:))];
        t_min    = [t_min min(t(:)) -max(t(:))];
        tperm(mask_P)    = tperm(mask_P) + 2*(t(mask_P) >= t0(mask_P));
        tperm(mask_N)    = tperm(mask_N) - 2*(t(mask_N) <= t0(mask_N));
        tfceperm(mask_P) = tfceperm(mask_P) + 2*(tfce(mask_P) >= tfce0(mask_P));
        tfceperm(mask_N) = tfceperm(mask_N) - 2*(tfce(mask_N) <= tfce0(mask_N));
      end
    else
      if n_cond == 1 % one-sample t-test
        label_matrix = [label_matrix; rand_label];
      else
        label_matrix = [label_matrix; rand_order_sorted];
      end

      if ~test_mode
        % maximum statistic
        tfce_max = [tfce_max max(tfce(:))];
        t_max    = [t_max max(t(:))];
        tfce_min = [tfce_min min(tfce(:))];
        t_min    = [t_min min(t(:))];
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
      stfce_min = sort(tfce_min,2,'descend');
      st_min    = sort(t_min,2,'descend');
  
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
        h1 = axes('position',[0 0 1 0.95],'Parent',Fgraph,'Visible','off');
        plot_distribution(stfce_max, tfce_max_th, 'tfce', alpha, col, 1, tfce0_max, tfce0_min);
        if ~show_permuted_designmatrix
          plot_distribution(st_max, t_max_th, 't-value', alpha, col, 2, t0_max, t0_min);
        end
      end
    
      if numel(job.conspec.n_perm) > 1
        if i > n_perm_break
          if isempty(find(tfce0_max > tfce_max_th(50:end,1)))
            fprintf('No FWE-corrected suprathreshold value after %d permutations found\n', n_perm_break);
            i = n_perm;
          end
        end  
      end
    end % test_mode

    
    if use_half_permutations
      if ~rem(i,progress_step) || ~rem(i+1,progress_step)
        if ~test_mode, cg_progress('Set',i,Fgraph); end
        drawnow
      end
    else
      if ~rem(i,progress_step)
        if ~test_mode, cg_progress('Set',i,Fgraph); end
        drawnow
      end
    end
  
    if use_half_permutations  
      i = i + 2;
    else
      i = i + 1;
    end
  
  end
  
  if ~test_mode, cg_progress('Clear',Fgraph); end
  
  try
    delete(hStopButton);
    spm_print;
  end
  
  if ~test_mode
    % get correct number of permutations in case that process was stopped
    n_perm = length(tfce_max);

    %---------------------------------------------------------------
    % corrected threshold based on permutation distribution
    %---------------------------------------------------------------
  
    % allow thresholds depending on # of permutations
    n_alpha = 3;
    if n_perm < 1000, n_alpha = 2; end
    if n_perm <  100, n_alpha = 1; end
      
    %---------------------------------------------------------------
    % save unpermuted TFCE map
    %---------------------------------------------------------------
    name = sprintf('TFCE_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d %s',Ic, str_permutation_method);
    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,tfce0);
    else
      spm_write_vol(Vt,tfce0);
    end
  
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
  
    if ~isempty(mask_P)
      nPtfcelog10(mask_P) = -log10(nPtfce(mask_P));
    end
    if ~isempty(mask_N)
      nPtfce(mask_N) = -nPtfce(mask_N);
      nPtfcelog10(mask_N) =  log10(nPtfce(mask_N));
    end
    
    nPtfce(mask_0) = NaN;
    nPtfcelog10(mask_0) = NaN;

    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,nPtfcelog10);
    else
      spm_write_vol(Vt,nPtfcelog10);
    end
  
    %---------------------------------------------------------------
    % save uncorrected p-values for T
    %---------------------------------------------------------------
    name = sprintf('%s_log_p_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d %s',xCon.STAT,Ic, str_permutation_method);
  
    nPtlog10 = zeros(size(t0));
  
    % estimate p-values
    nPt = tperm/n_perm;
   
    if ~isempty(mask_P)
      nPtlog10(mask_P) = -log10(nPt(mask_P));
    end
    if ~isempty(mask_N)
      nPt(mask_N) = -nPt(mask_N);
      nPtlog10(mask_N) =  log10(nPt(mask_N));
    end
    
    nPt(mask_0) = NaN;
    nPtlog10(mask_0) = NaN;

    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,nPtlog10);
    else
      spm_write_vol(Vt,nPtlog10);
    end
  
    %---------------------------------------------------------------
    % save corrected p-values for TFCE
    %---------------------------------------------------------------
    fprintf('Save corrected p-values.\n');

    name = sprintf('TFCE_log_pFWE_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d FWE %s',Ic, str_permutation_method);

    corrP = zeros(size(t));
  
    if ~isempty(mask_P)
      for t2 = tfce_max
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP(mask_P) = corrP(mask_P) + (t2 > tfce0(mask_P)  - tol);
      end
    end
    
    if ~isempty(mask_N)
      for t2 = tfce_min
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP(mask_N) = corrP(mask_N) - (t2 < tfce0(mask_N) + tol);
      end
    end
    
    corrP = corrP/n_perm;  
    corrPlog10 = zeros(size(tfce0));

    if ~isempty(mask_P)
      corrPlog10(mask_P) = -log10(corrP(mask_P));
    end
    
    if ~isempty(mask_N)
      corrP(mask_N) = -corrP(mask_N);
      corrPlog10(mask_N) =  log10(corrP(mask_N));
    end
    
    corrP(mask_0) = NaN;
    corrPlog10(mask_0) = NaN;

    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,corrPlog10);
    else
      spm_write_vol(Vt,corrPlog10);
    end
  
    %---------------------------------------------------------------
    % save corrected p-values for T
    %---------------------------------------------------------------
    name = sprintf('%s_log_pFWE_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d FWE %s',xCon.STAT,Ic, str_permutation_method);

    corrP = zeros(size(t));
  
    if ~isempty(mask_P)
      for t2 = t_max
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP(mask_P) = corrP(mask_P) + (t2 > t0(mask_P) - tol);
      end
    end
    if ~isempty(mask_N)
      for t2 = t_min
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP(mask_N) = corrP(mask_N) - (t2 < t0(mask_N) + tol);
      end
    end
    
    corrP = corrP / n_perm;  
    corrPlog10 = zeros(size(tfce0));

    if ~isempty(mask_P)
      corrP(mask_P) = corrP(mask_P);
      corrPlog10(mask_P) = -log10(corrP(mask_P));
    end
    if ~isempty(mask_N)
      corrP(mask_N) = -corrP(mask_N);
      corrPlog10(mask_N) =  log10(corrP(mask_N));
    end

    corrP(mask_0) = NaN;
    corrPlog10(mask_0) = NaN;
  
    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,corrPlog10);
    else
      spm_write_vol(Vt,corrPlog10);
    end
  
    %---------------------------------------------------------------
    % save corrected FDR-values for TFCE
    %---------------------------------------------------------------
    fprintf('Save corrected FDR-values.\n');
  
    name = sprintf('TFCE_log_pFDR_%04d',Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('TFCE %04d FDR %s',Ic, str_permutation_method);

    corrPfdr = NaN(size(t));
    corrPfdrlog10 = zeros(size(tfce0));

    if ~isempty(mask_P)
      [snP_pos,I_pos] = sort(nPtfce(mask_P));
      if ~isempty(snP_pos)
        corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
        corrPfdr_pos(I_pos) = corrPfdr_pos;
        corrPfdr(mask_P) = corrPfdr_pos;
        corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
      end
    end
  
    if ~isempty(mask_N)
      [snP_neg,I_neg] = sort(nPtfce(mask_N));
      if ~isempty(snP_neg)
        corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
        corrPfdr_neg(I_neg) = corrPfdr_neg;
        corrPfdr(mask_N) = corrPfdr_neg;
        corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
      end
    end
      
    corrPfdrlog10(mask_0) = NaN;

    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,corrPfdrlog10);
    else
      spm_write_vol(Vt,corrPfdrlog10);
    end
  
    %---------------------------------------------------------------
    % save corrected FDR-values for T
    %---------------------------------------------------------------
    name = sprintf('%s_log_pFDR_%04d',xCon.STAT,Ic);
    Vt.fname = fullfile(cwd,[name file_ext]);
    Vt.descrip = sprintf('%s %04d FDR %s',xCon.STAT,Ic, str_permutation_method);

    corrPfdr = NaN(size(t));
    corrPfdrlog10 = zeros(size(tfce0));

    if ~isempty(mask_P)
      if ~isempty(snP_pos)
        [snP_pos,I_pos] = sort(nPt(mask_P));
        corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
        corrPfdr_pos(I_pos) = corrPfdr_pos;
        corrPfdr(mask_P) = corrPfdr_pos;
        corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
      end
    end
  
    if ~isempty(mask_N)
      [snP_neg,I_neg] = sort(nPt(mask_N));
      if ~isempty(snP_neg)
        corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
        corrPfdr_neg(I_neg) = corrPfdr_neg;
        corrPfdr(mask_N) = corrPfdr_neg;
        corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
      end
    end
  
    corrPfdrlog10(mask_0) = NaN;

    if spm12
      Vt = spm_data_hdr_write(Vt);
      spm_data_write(Vt,corrPfdrlog10);
    else
      spm_write_vol(Vt,corrPfdrlog10);
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

  [hmax, xmax] = hist(val_max, 20);
    
  subplot(2,2,(2*order)-1)
  
  h = bar(xmax,hmax);
  set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);

  avg_h = mean(hmax);
  max_h = max(hmax);
  lim_x = xlim;

  % plot maximum observed value for unpermuted model
  hl = line([val0_max val0_max], [0 max_h]);
  set(hl,'Color',[0.3333 1 0],'LineWidth',2);
  text(0.95*lim_x(2),0.95*max_h,'Max. observed value ',...
    'Color',[0.3333 1 0],'HorizontalAlignment','Right','FontSize',8)
    
  % plot sign-flipped minimum observed value for unpermuted model
  if val0_min < 0
    hl = line([-val0_min -val0_min], [0 max_h]);
    set(hl,'Color',[0 0.6667 1],'LineWidth',2);
    text(0.95*lim_x(2),0.85*max_h,'Max. observed value (inverse contrast) ',...
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
    title(['Corr. threshold of ' name],'FontWeight','bold')
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
n_data = size(xX.X,1);

xKXs = spm_sp('Set',xX.W*xX.X);
xKXs.X = full(xKXs.X);
pKX = spm_sp('x-',xKXs);

Beta = Y*pKX';
res0 = Beta*single(xKXs.X') - Y;     %-Residuals
ResSS = double(sum(res0.^2,2));      %-Residual SSQ

trRV = n_data - rank(xX.X);
ResMS = ResSS/trRV;

T = zeros(dim);

if strcmp(xCon.STAT,'T')
  Bcov = pKX*pKX';
  con = Beta*c;

  T(ind_mask) = con./(eps+sqrt(ResMS*(c'*Bcov*c)));
else
  %-Compute ESS
  % Residual (in parameter space) forming matrix
  h  = spm_FcUtil('Hsqr',xCon,xKXs);
  
  ess = sum((h*Beta').^2,1)';
  MVM = ess/xCon.eidf;

  T(ind_mask) = MVM./ResMS;
end

%---------------------------------------------------------------

function xX = correct_xX(xX)

% vector of covariates and nuisance variables
iCG = [xX.iC xX.iG];
iHB = [xX.iH xX.iB];

% set columns with covariates and nuisance variables to zero
X = xX.X;
X(:,iCG) = 0;

ncol = size(X,2);

% calculate sum of columns
% The idea behind this is that for each factor the sum of all of its columns should be "1".
Xsum = zeros(size(X));
for i=1:ncol
  % only sum up columns without covariates and nuisance variables
  if isempty(find(iCG==i))
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
  
  % not sure whether this will always work but usually iB (subject effects) should be longer than iH (time effects)
  if length(iF{1}) > length(iF{2})
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
if isnan(df2(1)),
    
  if df1 == 0,
      
    % If df1 is zero, this is already a z-stat (this is here more for
    % compatibility).
    Z = G;
      
  elseif df1 == 1,
      
    % If rank(C) = 1, i.e., df1 = 1, this is r, so
    % do a Fisher's r-to-z stransformation
    Z = atanh(G);
      
  elseif df1 > 1,
      
    % If rank(C) > 1, i.e., df1 > 1, this is R^2, so
    % use a probit transformation.
    Z = -erfcinv(2*G)*sqrt(2); %Z = norminv(G);
      
  end

else
  siz = size(G);
  Z = zeros(siz);
  df2 = bsxfun(@times,ones(siz),df2);
  if df1 == 1,
    
    % Deal with precision issues working on each
    % tail separately
    idx = G > 0;
    %Z( idx) = -erfinv(2*palm_gcdf(-G( idx),1,df2( idx))-1)*sqrt(2);
    %Z(~idx) =  erfinv(2*palm_gcdf( G(~idx),1,df2(~idx))-1)*sqrt(2);
    Z( idx) =  erfcinv(2*palm_gcdf(-G( idx),1,df2( idx)))*sqrt(2);
    Z(~idx) = -erfcinv(2*palm_gcdf( G(~idx),1,df2(~idx)))*sqrt(2);
    
  elseif df1 == 0,
    
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

if df1 > 1,
  
  % G or F
  df2 = bsxfun(@times,ones(size(G)),df2);
  B = (df1.*G./df2)./(1+df1.*G./df2);
  gcdf = betainc(B,df1/2,df2/2);

elseif df1 == 1,
  
  % Student's t, Aspin's v
  df2 = bsxfun(@times,ones(size(G)),df2);
  ic = df2 == 1;
  in = df2 > 1e7;
  ig = ~(ic|in);
  gcdf = zeros(size(G));
  if any(ig(:)),
    gcdf(ig) = betainc(1./(1+G(ig).^2./df2(ig)),df2(ig)/2,.5)/2;
  end
  ig = G > 0 & ig;
  gcdf(ig) = 1 - gcdf(ig);
  if any(ic(:)),
    gcdf(ic) = .5 + atan(G(ic))/pi;
  end
  if any(in(:)),
    gcdf(ic) = palm_gcdf(G(in),0);
  end

elseif df1 == 0,
  
  % Normal distribution
  gcdf = erfc(-G/sqrt(2))/2;
  
elseif df1 < 0,
  
  % Chi^2, via lower Gamma incomplete for precision and speed
  %df2 = bsxfun(@times,ones(size(G)),df2);
  gcdf = palm_gammainc(G/2,df2/2,'lower');
  
end

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

if (isa (x, 'single') || isa (a, 'single') || isa (b, 'single'));
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

function eq = size_equal(a,b)
eq = isequal(size(a),size(b));

function a = iscomplex(X)
a = ~isreal(X);

function varargout = lgamma(varargin)
varargout{1:nargout} = gammaln(varargin{:});
