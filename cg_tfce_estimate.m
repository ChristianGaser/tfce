function cg_tfce_estimate(SPM, Ic, xCon)

% define stepsize for tfce
n_steps_tfce = 100;

% define histogram bins
n_hist_bins = 1000;

% colors and alpha levels
col = str2mat('b','g','r');
alpha = [0.05 0.01 0.001];

% Initialise random number generator
rand('state',sum(100*clock));

% load SPM file
if nargin < 1
  Pmat = spm_select(1,'SPM.mat','Select SPM.mat');
  load(Pmat)
  cwd = fileparts(Pmat);
else
  cwd = SPM.swd;
end
% analysis directory

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
if length(SPM.xX) > 1
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
  label = ones(1,n_subj);
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

% estimate # of permutations
% Anova/correlation: n_perm = (n1+n2+...+nk)!/(n1!*n2!*...*nk!)
if n_cond ~=1  % Anova/correlation
  n_perm = factorial(n_subj);
  for i=1:n_cond
    n_perm = n_perm/factorial(length(find(label == i)));
  end
else  % one-sample t-test: n_perm = 2^n
  n_perm = 2^n_subj;
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
n_perm = spm_input('How many permutations? ',1,'r',n_perm,1,[10 n_perm]);
vFWHM  = spm_input('Variance smoothing in FWHM (for low DFs) ','+1','e',0);

% compute unpermuted t-map
time0 = clock;

t0 = calc_glm(VY,X,c,Vmask,vFWHM,TH,W);

% calculate tfce of unpermuted t-map
tfce0 = tfceMex(t0, n_steps_tfce);
fprintf('Estimated time to run all permutations: %3.1f min\n',n_perm*etime(clock, time0)/60);

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

% check that label has correct dimension
sz = size(label);
if sz(1)>sz(2)
  label = label';
end

stopStatus = false;
i = 0;
spm_progress_bar('Init',n_perm,'Calculating','Permutations')

% only update progress bar every 1%
progress_step = max([1 round(n_perm/100)]);

tic
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
      rand_label = ones(1,n_subj);
    else % correlation or Anova
      rand_order = 1:n_subj;
      rand_label = label(rand_order);
    end
  else
    % init permutation and
    % check that each permutation is used only once
    if n_cond == 1 % one-sample t-test
      rand_label = sign(randn(1,n_subj));
      while any(ismember(rand_vector,rand_label,'rows'))
        rand_label = sign(randn(1,n_subj));
      end
    else % correlation or Anova
      rand_order = randperm(n_subj);
      rand_label = label(rand_order);
      while any(ismember(rand_vector,rand_label,'rows'))
        rand_order = randperm(n_subj);
        rand_label = label(rand_order);
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
    for j=1:n_subj
      Xperm(j,ind_con) = rand_label(j)*Xperm(j,ind_con);
    end
  else % correlation or Anova
    Xperm(:,ind_con) = Xperm(rand_order,ind_con);
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
  tfce = tfceMex(t, n_steps_tfce);

  tfce_max = [tfce_max max(tfce(:))];
  t_max   = [t_max max(t(:))];
  
  tfce_hist = tfce_hist + hist(tfce(find(tfce>0)), tfce_bins);
  t_hist = t_hist + hist(t(find(t>0)), t_bins);
  
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

  plot_distribution(tfce_max, tfce_max_th, 'tfce', alpha, col, 1, tfce0_min, tfce0_max);
  plot_distribution(t_max ,t_max_th, 't-value', alpha, col, 2, t0_min, t0_max);

  i = i + 1;

  drawnow
  
  if ~rem(i,progress_step)
    spm_progress_bar('Set',i);
  end
  
end
toc

%save rand_vector rand_vector
%save tfce_max tfce_max

spm_progress_bar('Clear')

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

% save TFCE map

name = sprintf('spmTFCE_%04d',Ic);
Vt.fname = fullfile(cwd,[name '.img']);
if vFWHM > 0
  Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d.img',vFWHM,Ic);
else
  Vt.descrip = sprintf('TFCE Contrast %04d.img',Ic);
end
spm_write_vol(Vt,tfce0);

% save corrected p-values for TFCE
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
uncorrP = zeros(size(tfce0));

% estimate p-values only with 100 steps to save time
tfce_cumsum = cumsum(tfce_hist);
for j=n_hist_bins:-1:1
  ind = find(uncorrP==0);
  tmp = min(find(tfce_cumsum>=ceil(j/n_hist_bins*sum(tfce_hist))));
  indp = find(tfce0(ind) >= tfce_bins(tmp));
  uncorrP(ind(indp)) = j/n_hist_bins;
  j/n_hist_bins
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

% estimate p-values only with 100 steps to save time
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

for j=1:n_alpha

  fprintf('Final corrected tfce threshold = %6d (p<%.3f)\n',round(tfce_max_th_final(j)),alpha(j));
  fprintf('Final corrected voxel height threshold = %6.2f (p<%.3f)\n',t_max_th_final(j),alpha(j));

  if tfce_max_th_final(j) > max([-tfce0_min tfce0_max])
    fprintf('No corrected suprathreshold tfce found.\n');
  else
    Vt = VY(1);
    Vt.dt = [spm_type('float32') spm_platform('bigend')];
    Vt.pinfo(1:2) = [1 0];

    % threshold unpermuted positive tfce-values
    if tfce0_max > tfce_max_th_final(j)
     tfce_corrected = tfce0;
     tfce_corrected(find(tfce_corrected <= tfce_max_th_final(j))) = 0;

     name = sprintf('spmT_%04d_tfce%d_pcorr%.3f',Ic,round(tfce_max_th_final(j)),alpha(j));
     Vt.fname = fullfile(cwd,[name '.img']);
     if vFWHM > 0
       Vt.descrip = sprintf('tfce=%d; vFWHM=%.1fmm, Contrast %04d.img',...
         round(tfce_max_th_final(j)),vFWHM,Ic);
     else
       Vt.descrip = sprintf('tfce=%d; Contrast %04d.img',...
         round(tfce_max_th_final(j)),Ic);
     end
     spm_write_vol(Vt,tfce_corrected);
    end
    
    % threshold unpermuted negative tfce values
    if tfce0_min < -tfce_max_th_final(j)
     tfce_corrected = tfce0;
     tfce_corrected(find(tfce_corrected >= -tfce_max_th_final(j))) = 0;

     name = sprintf('spmT_%04d_inverse_tfce%d_pcorr%.3f',Ic,round(tfce_max_th_final(j)),alpha(j));
     Vt.fname = fullfile(cwd,[name '.img']);
     if vFWHM > 0
       Vt.descrip = sprintf('tfce=%d; vFWHM=%.1fmm, inverse contrast %04d.img',...
         round(tfce_max_th_final(j)),vFWHM,Ic);
     else
       Vt.descrip = sprintf('tfce=%d; inverse contrast %04d.img',...
         round(tfce_max_th_final(j)),Ic);
     end
     spm_write_vol(Vt,tfce_corrected);
    end

  end

  fprintf('\n')
end

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

  % plot minimum observed value for unpermuted model
if 0
  hl = line([-val0_min -val0_min], [0 max_h]);
  set(hl,'Color','c','LineWidth',2);
  text(1.05*lim_x(1),0.9*max_h,'Min. observed negative value ',...
    'Color','c','HorizontalAlignment','Left','FontSize',8)
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
