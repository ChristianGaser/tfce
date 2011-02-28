function cg_tfce_estimate(SPM, Ic, xCon, n_perm_max, vFWHM, n_perm_break)

% define stepsize for tfce
n_steps_tfce = 100;

% colors and alpha levels
col   = [0 0 1; 0 0.5 0; 1 0 0];
alpha = [0.05 0.01 0.001];

% give same results each time
rand('state',0);

% tolerance for comparison
tol = 1e-4;	% Tolerance for comparing real numbers

% load SPM file
if nargin < 1
  Pmat = spm_select(1,'SPM.mat','Select SPM.mat');
  load(Pmat)
  cwd = fileparts(Pmat);
else
  cwd = SPM.swd;
end

%-Check that model has been estimated
try
  xX  = SPM.xX;         %-Design definition structure
  XYZ = SPM.xVol.XYZ;      %-XYZ coordinates
  S   = SPM.xVol.S;       %-search Volume {voxels}
  R   = SPM.xVol.R;       %-search Volume {resels}
  M   = SPM.xVol.M(1:3,1:3);   %-voxels to mm matrix
  VOX = sqrt(diag(M'*M))';    %-voxel dimensions
catch
  str = { 'This model has not been estimated.';...
          'Would you like to estimate it now?'};
  if spm_input(str,1,'bd','yes|no',[1,0],1)
     SPM = spm_spm(SPM);
     xX  = SPM.xX;          %-Design definition structure
     XYZ = SPM.xVol.XYZ;       %-XYZ coordinates
     S   = SPM.xVol.S;        %-search Volume {voxels}
     R   = SPM.xVol.R;        %-search Volume {resels}
     M   = SPM.xVol.M(1:3,1:3);   %-voxels to mm matrix
     VOX = sqrt(diag(M'*M))';    %-voxel dimensions
  else
    return
  end
end

% check that no temporal filter was used
if isstruct(SPM.xX.K)
  error('No first level analysis with temporal correlations allowed');
end

% whitening matrix
if isfield(SPM.xX,'W')
  W = SPM.xX.W;
  if any(W>1)
    warning('Whitening of the data is not yet supported');
    W = speye(size(SPM.xX.X,1));
  end
else
  W = speye(size(SPM.xX.X,1));
end

% design matrix
xX = SPM.xX;
X  = SPM.xX.X;

% threshold vector
TH = SPM.xM.TH;

if nargin < 3
  [Ic,xCon] = spm_conman(SPM,'T',1,...
        '  Select contrasts...',' for conjunction',1);
end

if length(Ic) > 1
  error('No conjunction allowed.');
end

% get contrast and find zero values in contrast
c       = xCon(Ic).c;
ind_con = find(c~=0);
c_name  = deblank(xCon(Ic).name);

% find exchangeability blocks using contrasts without zero values
unique_con   = unique(c(ind_con));
n_unique_con = length(unique_con);

% check for exchangeability blocks and design matrix
% maximal two exchangeability blocks allowed
switch n_unique_con
  case 1 % check whether the contrast is defined at columns for condition effects
    n_cond = length(find(SPM.xX.iH==ind_con));
    use_half_permutations = 1;
  case 2 % check whether the contrast is defined at columns for condition effects
    n_cond = 0;
    n_subj_cond = [];
    for j=1:n_unique_con
      c_unique_con = find(c==unique_con(j));
      for k=1:length(c_unique_con)
        n_cond = n_cond + length(find(SPM.xX.iH==c_unique_con(k)));
        n_subj_cond = [n_subj_cond sum(SPM.xX.X(:,find(SPM.xX.iH==c_unique_con(k))))];
      end
    end
    
    % check if sample size is equal for both conditions
    if n_subj_cond(1) == n_subj_cond(2)
      use_half_permutations = 1;
      disp('Equal sample sizes: half of permutations could be used.');
    else
      use_half_permutations = 0;
    end
    
    % comparison of two regressors (=interaction) not fully tested
    if n_cond == 0
      warning('Interaction between two regressors should work, but is not yet fully tested.')
    end
  otherwise
    error('Maximal two exchangeability blocks allowed.')
end

ind_unique_con = cell(n_unique_con,1);
for j=1:n_unique_con
  ind_unique_con{j} = find(c==unique_con(j));
end

n_subj = size(X,1);

% check design
switch n_cond
case 0 % correlation
  fprintf('Correlational design found\n');
  label = 1:n_subj;
case 1 % one-sample t-test
  fprintf('One sample t-test found\n');
  % use exchangeability blocks for labels
  label = zeros(1,n_subj);
  for j=1:n_unique_con
    for k=1:length(ind_unique_con{j})
      label(find(xX.X(:,ind_unique_con{j}(k))==1)) = j;
    end
  end
otherwise  % Anova with at least 2 groups
  fprintf('Anova found\n');
  % use exchangeability blocks for labels
  label = zeros(1,n_subj);
  for j=1:n_unique_con
    for k=1:length(ind_unique_con{j})
      label(find(xX.X(:,ind_unique_con{j}(k))==1)) = j;
    end
  end
end

fprintf('\n')

% get index for label values > 0
ind_label = find(label > 0);

n_subj_with_contrast = length(find(label > 0));

% estimate # of permutations
% Anova/correlation: n_perm = (n1+n2+...+nk)!/(n1!*n2!*...*nk!)
if n_cond ~=1  % Anova/correlation
  n_perm_full = factorial(n_subj_with_contrast);
  for i=1:n_cond
    n_perm_full = n_perm_full/factorial(length(find(label == i)));
  end
  n_perm_full = round(n_perm_full);
else  % one-sample t-test: n_perm = 2^n
  n_perm_full = 2^n_subj_with_contrast;
end

VY = SPM.xY.VY;

% load mask file
maskname = fullfile(cwd,'mask.img');
if ~exist(maskname)
  maskname = spm_select(1,'image','select mask image');
end
Vmask = spm_vol(maskname);

% if first image was not found you have to select all files again
if ~exist(VY(1).fname);
  P = spm_select(size(SPM.xY.VY,1),'image','select images');
  VY = spm_vol(P);
  SPM.xY.VY = VY;
  
  % update SPM
  save(Pmat,'SPM');
end
clear SPM

% choose number of permutations
if nargin < 4
  n_perm = spm_input('How many permutations? ',1,'r',n_perm_full,1,[10 n_perm_full]);
else
  n_perm = min([n_perm n_perm_max]);
end
if nargin < 5
  vFWHM  = spm_input('Variance smoothing (for low DFs) ','+1','e',0);
end

W = sparse(eye(length(W)));
warning('Whitening is not considered!');

% compute unpermuted t-map
t0 = calc_glm(VY,X,c,Vmask,vFWHM,TH,W);

% calculate tfce of unpermuted t-map
try
    tfce0 = tfceMex(t0, n_steps_tfce,1);
catch
    tfce0 = tfceMex_noopenmp(t0, n_steps_tfce);
end

% get largest tfce
tfce0_min = min(tfce0(:));
tfce0_max = max(tfce0(:));
t0_min  = min(t0(:));
t0_max  = max(t0(:));

% prepare countings
t_max        = [];
t_max_th     = [];
t_th         = [];
tfce_max     = [];
tfce_max_th  = [];
tfce_th      = [];
label_matrix = [];

% general initialization
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
figure(Fgraph)

h = axes('position',[0.45 0.95 0.1 0.05],'Units','normalized','Parent',...
    Fgraph,'Visible','off');
text(0.5,0.5,c_name,...
    'FontSize',spm('FontSize',16),...
    'FontWeight','Bold',...
    'FontName',spm_platform('Font','times'),...
    'HorizontalAlignment','Center',...
    'VerticalAlignment','middle')

% check that label has correct dimension
sz = size(label);
if sz(1)>sz(2)
  label = label';
end

stopStatus = false;
spm_progress_bar('Init',n_perm,'Calculating','Permutations')
cg_progress('Init',n_perm,'Calculating','Permutations')

% update interval for progress bar
progress_step = max([1 round(n_perm/1000)]);

for i=1:n_perm

  % randomize subject vector
  if i==1 % first permutation is always unpermuted model
    if n_cond == 1 % one-sample t-test
      rand_label = ones(1,n_subj_with_contrast);
    else % correlation or Anova
      rand_order = 1:n_subj_with_contrast;
      rand_label = label(ind_label(rand_order));
    end
  else
    % init permutation and
    % check that each permutation is used only once
    if n_cond == 1 % one-sample t-test
      rand_label = sign(randn(1,n_subj_with_contrast));
      while any(ismember(label_matrix,rand_label,'rows'))
        rand_label = sign(randn(1,n_subj_with_contrast));
      end
    else % correlation or Anova
      rand_order = randperm(n_subj_with_contrast);
      rand_label = label(ind_label(rand_order));
      while any(ismember(label_matrix,rand_label,'rows'))
        rand_order = randperm(n_subj_with_contrast);
        rand_label = label(ind_label(rand_order));
      end
    end    
  end   
  
  % update label_matrix and order_matrix for checking of unique permutations
  label_matrix = [label_matrix; rand_label];

  % add Stop button after 20 iterations
  if i==20
    hStopButton = uicontrol(Fgraph,...
      'position',[10 10 70 20],...
      'style','toggle',...
      'string','Stop',...
      'backgroundcolor',[1 .5 .5]); % light-red
  end
  if i>=20
    stopStatus = get(hStopButton,'value');
  end
  
  % check Stop status
  if (stopStatus == true)
    fprintf('Stopped after %d iterations.\n',i);
    break; % stop the permutation loop
  end
    
  % change design matrix according to permutation order
  % only permute columns, where contrast is defined
  Xperm = X;
  
  if n_cond==1 % one-sample t-test
    % only change sign in the design matrix
    for j=1:n_subj_with_contrast
      Xperm(ind_label(j),ind_con) = rand_label(j)*Xperm(ind_label(j),ind_con);
    end
  else % correlation or Anova
    Xperm(ind_label,ind_con) = Xperm(ind_label(rand_order),ind_con);
  end

  % calculate permuted t-map
  if i==1
    t = t0;
    tfce = tfce0;
    nPt = ones(size(t));
    nPtfce = ones(size(t));
  else
    % -----------------------------------------------------
    % -----------------------------------------------------
    % What about W for whitening the data??? Should this also permuted???
    % -----------------------------------------------------
    % -----------------------------------------------------
    
    t = calc_glm(VY,Xperm,c,Vmask,vFWHM,TH,W);
    
    % compute tfce
    try
      tfce = tfceMex(t, n_steps_tfce);
    catch
      tfce = tfceMex_noopenmp(t, n_steps_tfce);
    end  
    
    % uncorrected p-values
    nPt = nPt + (t>=t0);
    nPtfce = nPtfce + (tfce>=tfce0);
  end


  tfce_max = [tfce_max max(tfce(:))];
  t_max   = [t_max max(t(:))];
    
  % use cummulated sum to find threshold
  tfce_max = sort(tfce_max);
  t_max    = sort(t_max);

  % find corrected thresholds
  ind_max  = ceil((1-alpha).*length(t_max));
  t_max_th = [t_max_th; t_max(ind_max);];
  
  ind_max     = ceil((1-alpha).*length(tfce_max));
  tfce_max_th = [tfce_max_th; tfce_max(ind_max);];
        
  % plot thresholds and histograms
  figure(Fgraph)
  h1 = axes('position',[0 0 1 0.95],'Parent',Fgraph,'Visible','off');
  plot_distribution(tfce_max, tfce_max_th, 'tfce', alpha, col, 1, tfce0_min, tfce0_max);
  plot_distribution(t_max ,t_max_th, 't-value', alpha, col, 2, t0_min, t0_max);

  drawnow

  % Check for suprathreshold values
  if nargin > 5
    if i > n_perm_break
      if isempty(find(tfce0_max > tfce_max_th(50:end,1)))
        fprintf('No FWE-corrected suprathreshold value after %d permutations found\n', n_perm_break);
        i = n_perm;
      end
    end
  end
    
  if ~rem(i,progress_step)
    cg_progress('Set',i)
    spm_progress_bar('Set',i);
  end
  
end

cg_progress('Clear')
spm_progress_bar('Clear')

spm_print

% get correct number of permutations in case that process was stopped
n_perm = length(tfce_max);

nPt = nPt/n_perm;
nPtfce = nPtfce/n_perm;

%---------------------------------------------------------------
% corrected threshold based on permutation distribution
%---------------------------------------------------------------

tfce_max_th_final = tfce_max_th(end,:);
t_max_th_final    = t_max_th(end,:);

% allow thresholds depending on # of permutations
n_alpha = 3;
if n_perm < 1000, n_alpha = 2; end
if n_perm <  100, n_alpha = 1; end

% pepare output files
Vt = VY(1);
Vt.dt(1)    = 16;
Vt.pinfo(1) = 1;

% create mask
mask_P = find(nPt<1);

% save unpermuted TFCE map
name = sprintf('spmTFCE_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,tfce0);

% save ascii file with number of permutations
fid = fopen(fullfile(cwd,[name '.txt']),'w');
fprintf(fid,'%d\n',n_perm);
fclose(fid);

% save corrected p-values for TFCE
fprintf('Save corrected p-values.\n');
corrP = ones(size(t));

for t2 = tfce_max
	%-FEW-corrected p is proportion of randomisation greater or
	% equal to statistic.
	%-Use a > b -tol rather than a >= b to avoid comparing
	% two reals for equality.
	corrP = corrP + (t2 > tfce0-tol);
end
corrP_vol = NaN(size(t));
corrP_vol(mask_P) = corrP(mask_P) / n_perm;  

name = sprintf('TFCE_log_pFWE_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,-log10(corrP_vol));

% save corrected p-values for T
corrP = ones(size(t));

for t2 = t_max
	%-FEW-corrected p is proportion of randomisation greater or
	% equal to statistic.
	%-Use a > b -tol rather than a >= b to avoid comparing
	% two reals for equality.
	corrP = corrP + (t2 > t0-tol);
end
corrP_vol = NaN(size(t));
corrP_vol(mask_P) = corrP(mask_P) / n_perm;  

name = sprintf('T_log_pFWE_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('T Contrast %04d.img',Ic);
end
spm_write_vol(Vt,-log10(corrP_vol));

% save uncorrected p-values for TFCE
fprintf('Save uncorrected p-values.\n');

name = sprintf('TFCE_log_p_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,-log10(nPtfce));

% save uncorrected p-values for T
name = sprintf('T_log_p_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('T Contrast %04d.img',Ic);
end
spm_write_vol(Vt,-log10(nPt));

% save corrected FDR-values for TFCE
fprintf('Save corrected FDR-values.\n');

[snP_pos,I_pos]=sort(nPtfce(mask_P));
corrPfdr_pos=snpm_P_FDR([],[],'P',[],snP_pos);
corrPfdr_pos(I_pos) = corrPfdr_pos;
corrPfdr_pos_vol = NaN(size(t));
corrPfdr_pos_vol(mask_P) = corrPfdr_pos;

name = sprintf('TFCE_log_pFDR_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,-log10(corrPfdr_pos_vol));

% save corrected FDR-values for T
[snP_pos,I_pos]=sort(nPt(mask_P));
corrPfdr_pos=snpm_P_FDR([],[],'P',[],snP_pos);
corrPfdr_pos(I_pos) = corrPfdr_pos;
corrPfdr_pos_vol = NaN(size(t));
corrPfdr_pos_vol(mask_P) = corrPfdr_pos;

name = sprintf('T_log_pFDR_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('T Contrast %04d.img',Ic);
end
spm_write_vol(Vt,-log10(corrPfdr_pos_vol));

return
%---------------------------------------------------------------

function plot_distribution(val_max,val_th,name,alpha,col,order,val0_min,val0_max)

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

  [hmax, xmax] = hist(val_max, 10);
    
  subplot(2,2,(2*order)-1)
  
  bar(xmax,hmax)

  avg_h = mean(hmax);
  max_h = max(hmax);
  lim_x = xlim;

  % plot maximum observed value for unpermuted model
  hl = line([val0_max val0_max], [0 max_h]);
  set(hl,'Color','m','LineWidth',2);
  text(0.95*lim_x(2),0.95*max_h,'Max. observed value ',...
    'Color','m','HorizontalAlignment','Right','FontSize',8)
    
  % plot thresholds
  for j=1:n_alpha
    hl = line([val_th(n,j) val_th(n,j)], [0 max_h]);
    set(hl,'Color',col(j,:),'LineStyle','--');
    text(0.95*lim_x(2),(0.5+0.1*j)*max_h,['p<' num2str(alpha(j))],...
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
    semilogy(1:n,val_th(1:n,:))
    yl = log10(ylim);
    ylim(10.^[floor(yl(1)) ceil(yl(2))])
  else
    plot(1:n,val_th(1:n,:))
  end
  if corr
    title(['Corr. threshold of ' name],'FontWeight','bold')
  else
    title(['Uncorr. threshold of ' name],'FontWeight','bold')
  end
  ylabel('Threshold')
  xlabel('Permutations')   
end
return

%---------------------------------------------------------------
function t = calc_glm(V,X,c,Vmask,vFWHM,TH,W)
% compute t-statistic using GLM
%
% V   - memory mapped files
% X   - design structure
% c   - contrast
% Vmask - memory mapped mask image
% vFWHM - filter width for variance smoothing
% TH  - threshold vector
% W   - whitening matrix

n_subj = size(X,1);
n_beta = size(X,2);

X = W*X;
pKX  = pinv(X);
Bcov = pinv(X'*X);
trRV = n_subj - rank(X);

vx  = sqrt(sum(V(1).mat(1:3,1:3).^2));

dim = V(1).dim(1:3);

% try mex-file, which is faster
try
 [Beta, ResMS] = cg_glm_get_Beta_ResSS(V,Vmask,X,pKX,TH,W);
 ResMS = reshape(ResMS,dim);
catch
 [Beta, ResMS] = calc_beta(V,Vmask,X,pKX,TH,W);
end

Beta = reshape(Beta,[prod(dim) n_beta]);
ResMS = ResMS/trRV;

con = Beta*c;
clear Beta
con = reshape(con,dim);

% Code for variance smoothing was shameless taken from snpm3 code...
% Blurred mask is used to truncate kernal to brain; if not
% used variance at edges would be underestimated due to
% convolution with zero activity out side the brain.
if vFWHM > 0
  mask = spm_read_vols(Vmask);
  Q = find(mask > 0);
  SmResMS   = zeros(dim);
  SmMask    = zeros(dim);
  TmpVol    = zeros(dim);
  TmpVol(Q) = ones(size(Q));
  spm_smooth(TmpVol,SmMask,vFWHM./vx);
  TmpVol(Q) = ResMS(Q);
  spm_smooth(TmpVol,SmResMS,vFWHM./vx);
  ResMS(Q)  = SmResMS(Q)./SmMask(Q);
end

t = con./(eps+sqrt(ResMS*(c'*Bcov*c)));

return

%---------------------------------------------------------------

function [beta, ResSS] = calc_beta(VY,Vm,X,pKX,TH,W);

n = size(VY,1);
m = size(pKX,1);
beta = zeros([Vm.dim(1:3) m]);
ResSS = zeros(Vm.dim(1:3));

for j=1:Vm.dim(3),

 M  = spm_matrix([0 0 j]);
 mask = spm_slice_vol(Vm,M,Vm.dim(1:2),[1 0]);
 ind = find(mask>0);
 
 if ~isempty(ind)
  Y = zeros([length(ind) n],'single');
  pXY = zeros([prod(Vm.dim(1:2)) m],'single');
  res = zeros(Vm.dim(1:2));

  % Load slice j from all images
  for i=1:n
   tmp = spm_slice_vol(VY(i),M,Vm.dim(1:2),[1 0]);
   Y(:,i) = single(W(i,i)*tmp(ind));
  end
  
  pXY(ind,:) = Y*single(pKX'); 
  
  res0 = pXY(ind,:)*single(X') - Y; %-Residuals
  res(ind) = double(sum(res0.^2,2));
  ResSS(:,:,j) = res;   %-Residual SSQ

  pXY = reshape(pXY,[Vm.dim(1:2) m]);
  beta(:,:,j,:) = double(pXY);
 end
end

function [str] = estimate(t)
% Gives an appropriately unitized string of a duration in seconds.
%
% [STR] = ESTIMATE(T)
%

% License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
% 
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================
minutes = t./60;
hours = t./3600;
days = hours./24;

if days > 1
  str = sprintf('%d days %02.1f hr', floor(days),24*(days-floor(days)));
elseif hours > 1
  str = sprintf('%d:%02.0f hr', floor(hours),60*(hours-floor(hours)));
elseif minutes > 1
  str = sprintf('%d:%02.0f min',floor(minutes),60*(minutes-floor(minutes)));
else
  str = sprintf('%02.1f sec',t);
end
