function cg_tfce_estimate(SPM, Ic, xCon, n_perm_max, vFWHM, n_perm_break)

% define stepsize for tfce
n_steps_tfce = 100;

% define histogram bins
% use value > 1000 to reliable estimate p<0.001 levels
n_hist_bins = 1100;

% colors and alpha levels
col = str2mat('b','g','r');
alpha = [0.05 0.01 0.001];

% give same results each time
rand('state',0);

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
  XYZ  = SPM.xVol.XYZ;      %-XYZ coordinates
  S   = SPM.xVol.S;       %-search Volume {voxels}
  R   = SPM.xVol.R;       %-search Volume {resels}
  M   = SPM.xVol.M(1:3,1:3);   %-voxels to mm matrix
  VOX  = sqrt(diag(M'*M))';    %-voxel dimensions
catch
  str = { 'This model has not been estimated.';...
      'Would you like to estimate it now?'};
  if spm_input(str,1,'bd','yes|no',[1,0],1)
     SPM = spm_spm(SPM);
     xX  = SPM.xX;          %-Design definition structure
     XYZ = SPM.xVol.XYZ;       %-XYZ coordinates
     S  = SPM.xVol.S;        %-search Volume {voxels}
     R  = SPM.xVol.R;        %-search Volume {resels}
     M  = SPM.xVol.M(1:3,1:3);   %-voxels to mm matrix
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
    warning('Whitening of the data is probably not working');
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
c = xCon(Ic).c;
ind_con = find(c~=0);
c_name = deblank(xCon(Ic).name);

% find exchangeability blocks using contrasts without zero values
unique_con  = unique(c(ind_con));
n_unique_con = length(unique_con);

% check for exchangeability blocks and design matrix
% maximal two exchangeability blocks allowed
switch n_unique_con
  case 1 % check whether the contrast is defined at columns for condition effects
    n_cond = length(find(SPM.xX.iH==ind_con));
  case 2 % check whether the contrast is defined at columns for condition effects
    n_cond = 0;
    for j=1:n_unique_con
     c_unique_con = find(c==unique_con(j));
     for k=1:length(c_unique_con)
       n_cond = n_cond + length(find(SPM.xX.iH==c_unique_con(k)));
     end
    end
    % comparison of two regressors (=interaction) not fully tested
    if n_cond == 0
      warning('Interaction between two regressors was not yet fully tested.')
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
  n_perm = factorial(n_subj_with_contrast);
  for i=1:n_cond
    n_perm = n_perm/factorial(length(find(label == i)));
  end
  n_perm = round(n_perm);
else  % one-sample t-test: n_perm = 2^n
  n_perm = 2^n_subj_with_contrast;
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
end
clear SPM

% choose number of permutations
if nargin < 4
  n_perm = spm_input('How many permutations? ',1,'r',n_perm,1,[10 n_perm]);
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

% get vector for histogram bins
tfce_bins = linspace(0, max(abs(tfce0(:))), n_hist_bins);
t_bins    = linspace(0, max(abs(t0(:))), n_hist_bins);

% prepare countings
t_hist = zeros(1, n_hist_bins);
tfce_hist = zeros(1, n_hist_bins);
t_max    = [];
t_max_th  = [];
t_th = [];
tfce_max  = [];
tfce_max_th = [];
tfce_th = [];
rand_vector = [];

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
i = 0;
spm_progress_bar('Init',n_perm,'Calculating','Permutations')
cg_progress('Init',n_perm,'Calculating','Permutations')

% update interval for progress bar
progress_step = max([1 round(n_perm/1000)]);

while(i < n_perm)

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

  % randomize subject vector
  if i==0 % first permutation is always unpermuted model
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
      while any(ismember(rand_vector,rand_label,'rows'))
        rand_label = sign(randn(1,n_subj_with_contrast));
      end
    else % correlation or Anova
      rand_order = randperm(n_subj_with_contrast);
      rand_label = label(ind_label(rand_order));
      while any(ismember(rand_vector,rand_label,'rows'))
        rand_order = randperm(n_subj_with_contrast);
        rand_label = label(ind_label(rand_order));
      end
    end    
  end   
  
  % update rand_vector for checking of unique permutations
  rand_vector = [rand_vector; rand_label];
  
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
  if i==0
    t = t0;
  else
    % -----------------------------------------------------
    % -----------------------------------------------------
  % What about W for whitening the data??? Should this also permuted???
    % -----------------------------------------------------
    % -----------------------------------------------------
    t = calc_glm(VY,Xperm,c,Vmask,vFWHM,TH,W);
  end

  % compute tfce
  try
    tfce = tfceMex(t, n_steps_tfce);
  catch
    tfce = tfceMex_noopenmp(t, n_steps_tfce);
  end  

  tfce_max = [tfce_max max(tfce(:))];
  t_max   = [t_max max(t(:))];
  
  tfce_gt0 = tfce(find(tfce>0));
  if ~isempty(tfce_gt0)
      tfce_hist = tfce_hist + hist(tfce_gt0, tfce_bins);
  end
  t_gt0 = t(find(tfce>0));
  if ~isempty(t_gt0)
      t_hist = t_hist + hist(t(find(t>0)), t_bins);
  end
  
  % use cummulated sum to find threshold
  tfce_max = sort(tfce_max);
  t_max   = sort(t_max);

  % find corrected thresholds
  ind_max = ceil((1-alpha).*length(t_max));
  t_max_th = [t_max_th; t_max(ind_max);];
  
  ind_max = ceil((1-alpha).*length(tfce_max));
  tfce_max_th = [tfce_max_th; tfce_max(ind_max);];
        
  % plot thresholds and histograms
  figure(Fgraph)
  h1 = axes('position',[0 0 1 0.95],'Parent',Fgraph,'Visible','off');
  plot_distribution(tfce_max, tfce_max_th, 'tfce', alpha, col, 1, tfce0_min, tfce0_max);
  plot_distribution(t_max ,t_max_th, 't-value', alpha, col, 2, t0_min, t0_max);

  drawnow

  i = i + 1;

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

%---------------------------------------------------------------
% corrected threshold based on permutation distribution
%---------------------------------------------------------------

tfce_max_th_final = tfce_max_th(end,:);
t_max_th_final  = t_max_th(end,:);

% allow thresholds depending on # of permutations
n_alpha = 3;
sz_val_max = length(tfce_max);
if sz_val_max < 1000, n_alpha = 2; end
if sz_val_max <  100, n_alpha = 1; end

n_perm = length(tfce_max);

% pepare output files
Vt = VY(1);
Vt.dt(1) = 16;
Vt.pinfo(1) = 1;

% save unpermuted TFCE map
name = sprintf('spmTFCE_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,tfce0);
fid = fopen(fullfile(cwd,[name '.txt']),'w');
fprintf(fid,'%d\n',n_perm);
fclose(fid);

% save corrected p-values for TFCE
fprintf('Save corrected p-values.\n');
corrP = zeros(size(tfce0));

for j=n_perm:-1:1
  ind = find(corrP==0);
  indp = find(tfce0(ind) >= tfce_max(j));
  corrP(ind(indp)) = j/n_perm;
end

name = sprintf('TFCE_corrP_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,corrP);

% save corrected p-values for T
corrP = zeros(size(t0));

for j=n_perm:-1:1
  ind = find(corrP==0);
  indp = find(t0(ind) >= t_max(j));
  corrP(ind(indp)) = j/n_perm;
end

name = sprintf('T_corrP_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('T Contrast %04d.img',Ic);
end
spm_write_vol(Vt,corrP);

% save uncorrected p-values for TFCE
fprintf('Save uncorrected p-values.\n');
uncorrP = zeros(size(tfce0));

% estimate p-values
tfce_cumsum = cumsum(tfce_hist);
for j=n_hist_bins:-1:1
  ind = find(uncorrP==0);
  tmp = min(find(tfce_cumsum>=ceil(j/n_hist_bins*sum(tfce_hist))));
  indp = find(tfce0(ind) >= tfce_bins(tmp));
  uncorrP(ind(indp)) = j/n_hist_bins;
end

name = sprintf('TFCE_P_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,uncorrP);

% save uncorrected p-values for T
uncorrP = zeros(size(t0));

% estimate p-values
t_cumsum = cumsum(t_hist);
for j=n_hist_bins:-1:1
  ind = find(uncorrP==0);
  tmp = min(find(t_cumsum>=ceil(j/n_hist_bins*sum(t_hist))));
  indp = find(t0(ind) >= t_bins(tmp));
  uncorrP(ind(indp)) = j/n_hist_bins;
end

name = sprintf('T_P_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('T Contrast %04d.img',Ic);
end
spm_write_vol(Vt,uncorrP);

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
  alpha = alpha(1:n_alpha);
  val_th = val_th(:,1:n_alpha);

  [hmax, xmax] = hist(val_max, 10);
    
  subplot(2,2,(2*order)-1)
  
  bar(xmax,hmax)
  avg_h = mean(hmax);
  max_h = max(hmax);
  lim_x = xlim;
  for j=1:n_alpha
    hl = line([val_th(n,j) val_th(n,j)], [0 max_h]);
    set(hl,'Color',col(j));
    text(0.95*lim_x(2),(0.5+0.1*j)*max_h,['p<' num2str(alpha(j))],...
      'Color',col(j),'HorizontalAlignment','Right','FontSize',8)
  end

  % plot maximum observed value for unpermuted model
  hl = line([val0_max val0_max], [0 max_h]);
  set(hl,'Color','m','LineWidth',2);
  text(0.95*lim_x(2),0.95*max_h,'Max. observed value ',...
    'Color','m','HorizontalAlignment','Right','FontSize',8)
  
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
  SmResMS  = zeros(dim);
  SmMask  = zeros(dim);
  TmpVol  = zeros(dim);
  TmpVol(Q) = ones(size(Q));
  spm_smooth(TmpVol,SmMask,vFWHM./vx);
  TmpVol(Q) = ResMS(Q);
  spm_smooth(TmpVol,SmResMS,vFWHM./vx);
  ResMS(Q) = SmResMS(Q)./SmMask(Q);
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
