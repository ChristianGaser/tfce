% Main loop
for con = 1:length(Ic0)
  Ic = Ic0(con);
  xCon = SPM.xCon(Ic);
  
  [isValid, n_perm, n_perm_break] = initializeAndCheckErrors(job, Ic, SPM);
  if ~isValid
      fprintf('Initialization error.\n');
      return
  end

  [isValid, exch_blocks, n_exch_blocks, is_eoi, F_contrast_multiple_rows, use_half_permutations, c0] = calculateContrastsAndExchangeabilityBlocks(xCon, xX);
  if ~isValid
      fprintf('Initialization error.\n');
      return
  end

  fprintf('Use contrast #%d of %s\n',Ic,job.data{1})

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

  [Rz, str_permutation_method, nuisance_method] = initializeNuisanceMethod(nuisance_method, xX, c0, interaction_design);

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
      
      [t0, df2, SmMask] = calc_GLM_voxelwise(Y,xX,SPM.xC(voxel_covariate),xCon,ind_mask,VY(1).dim,C,[],ind_X,pinv_method);
    else
      [t0, df2, SmMask] = calc_GLM(Y,xX,xCon,ind_mask,VY(1).dim,vFWHM);
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
      % use bilateral filter of t-map to increase SNR, see LISA paper (Lohmann et al. 2018)
      if filter_bilateral
        var_t0 = var(t0(find(t0~=0 & ~isnan(t0) & ~isinf(t0))));
        t0 = double(cat_vol_bilateral(single(t0),2,2,2,2,var_t0));
      end
      
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
      [permutedResults, stopStatus] = executePermutation(perm, n_perm, xX, label, ind_label, job, SPM, n_data_with_contrast, n_cond, n_unique_labels, ind_label_gt0);
      if stopStatus
          break;
      end

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

  processAndSaveResults(permutedResults, SPM, xCon, cwd);
end

function [permutedResults, stopStatus] = executePermutation(perm, n_perm, xX, label, ind_label, job, SPM, n_data_with_contrast, n_cond, n_unique_labels, ind_label_gt0, test_mode);

  permutedResults = []; % This should be structured according to what you expect to collect from each permutation
  
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
          stopStatus = true;
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
        t = calc_GLM(Y*(Pset'*Rz),xXperm,xCon,ind_mask,VY(1).dim,vFWHM,SmMask);
      else
        if voxel_covariate
          t = calc_GLM_voxelwise(Y,xXperm,SPM.xC(voxel_covariate),xCon,ind_mask,VY(1).dim,C,Pset,ind_X,pinv_method);
        else
          t = calc_GLM(Y,xXperm,xCon,ind_mask,VY(1).dim,vFWHM,SmMask);
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
        if filter_bilateral
          t = double(cat_vol_bilateral(single(t),2,2,2,2,var_t0));
        end
        
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
      if (perm >= stop_if_now_FWEeffects_found) && (tfce0_max < 0.95*stfce_max(ind_max(alpha==0.05)) && -tfce0_min < 0.95*stfce_max(ind_max(alpha==0.05)))
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
        fprintf('\nCorrelation between between parametric and non-parametric statistic is cc=%g, which means that your design and optionally your nuisance paramters were correctly recognized.\n',cc(1,2));
      end
      check_validity = true;
    end

  end % test_mode

end

function [Rz, str_permutation_method, nuisance_method] = initializeNuisanceMethod(nuisance_method, xX, c0, interaction_design)

  [indi, indj] = find(c0~=0);
  ind_X = unique(indi)';

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

end

function [isValid, n_perm, n_perm_break] = initializeAndCheckErrors(job, Ic, SPM, xCon)

  % Assume initialization is valid initially
  isValid = true;
  n_perm = job.conspec.n_perm(1);
  if numel(job.conspec.n_perm) > 1
      n_perm_break = job.conspec.n_perm(2);
  else
      n_perm_break = inf; % No break condition if only one permutation count is specified
  end

  % Example error check
  if length(Ic) > 1
      fprintf('ERROR: No conjunction allowed.\n');
      isValid = false;
      return
  end

  % Insert additional initialization and error checking as needed
end

function [isValid, exch_blocks, n_exch_blocks, is_eoi, F_contrast_multiple_rows, use_half_permutations, c0] = calculateContrastsAndExchangeabilityBlocks(xCon, xX)

  % Assume initialization is valid initially
  isValid = true;
  repeated_anova = ~isempty(xX.iB);

  F_contrast_multiple_rows = 0; % Default value

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
      isValid = false;
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

end
