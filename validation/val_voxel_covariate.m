function val_voxel_covariate
% Calibration and cost of the voxel-/vertex-wise covariate path.
%
% A covariate given as one image per subject makes the design matrix vary across
% the brain: at every element the covariate columns are replaced by the local
% values and the model is refitted. Where the contrast is defined on the
% voxel-wise covariate, the covariate images are permuted together with the
% design.
%
% Two scenarios are checked under the null:
%   (a) the voxel-wise covariate IS the effect of interest (it gets permuted)
%   (b) the voxel-wise covariate is a NUISANCE variable (it does not get permuted)
%
% The real calc_GLM_voxelwise of the toolbox is used, not a reimplementation.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

% The vectorised solve inside calc_GLM_voxelwise is also checked against the
% per-voxel pseudoinverse loop it replaces, which is kept in the toolbox as
% calc_GLM_voxelwise_loop and is used here as the reference.
%
val_util('header','Voxel-/vertex-wise covariates');
val_util('extract','calc_GLM_voxelwise','calc_GLM_voxelwise_batch', ...
         'calc_GLM_voxelwise_loop','calc_GLM','pinv2','svdecon');

rng(3);
n_vox = 3000; n_perm = 200;
dim = [n_vox 1 1]; ind_mask = (1:n_vox)';

fprintf('\n  uncorrected false-positive rate at alpha = 0.05 (nominal 0.050)\n');
fprintf('  %-44s %-10s\n','scenario','FPR');
fprintf('  %s\n', repmat('-',1,60));

for scenario = 1:3
  n = 30;

  % voxel-wise covariate: one positive value per subject and voxel
  C = 0.5 + 0.1*randn(n_vox, n);

  switch scenario
  case 1
    % (a) the covariate is the effect of interest
    name  = 'covariate of interest (permuted), t';
    X     = [ones(n,1) mean(C,1)'];      % column 2 = subject mean of the covariate
    c     = [0; 1];
    STAT  = 'T'; eidf = 1;
    ind_X = 2;
    cols  = 2;
  case 2
    % (b) the covariate is a nuisance variable, the group effect is tested
    name  = 'covariate as nuisance (not permuted), t';
    g     = [ones(15,1); zeros(15,1)];
    X     = [g 1-g mean(C,1)'];
    c     = [1; -1; 0];
    STAT  = 'T'; eidf = 1;
    ind_X = 1:2;
    cols  = 3;
  case 3
    % (c) F-test on a 3-group effect, covariate as nuisance
    name  = 'covariate as nuisance (not permuted), F';
    g     = repmat((1:3)', 10, 1);
    X     = double([g==1 g==2 g==3 mean(C,1)']);
    c     = [1 0; -1 1; 0 -1; 0 0];
    STAT  = 'F'; eidf = rank(c);
    ind_X = 1:3;
    cols  = 4;
  end

  xX.X = X; xX.W = speye(n);
  xC.cols = cols;
  xCon.c = c; xCon.STAT = STAT; xCon.eidf = eidf;

  % null data: independent of both the covariate and the group
  Y = randn(n_vox, n);

  t0 = calc_GLM_voxelwise(Y, xX, xC, xCon, ind_mask, dim, C, [], ind_X, 1);
  t0 = t0(ind_mask);

  cnt = zeros(n_vox,1);
  tic
  for p = 1:n_perm
    pr   = randperm(n);
    Pset = sparse(1:n, pr, 1, n, n);

    xXperm   = xX;
    xXperm.X = X;
    xXperm.X(:,ind_X) = Pset*X(:,ind_X);   % Draper-Stoneman, as in the toolbox

    tp = calc_GLM_voxelwise(Y, xXperm, xC, xCon, ind_mask, dim, C, Pset, ind_X, 1);
    cnt = cnt + (tp(ind_mask) >= t0);
  end
  el = toc;

  nPt = cnt/n_perm;
  fpr = mean(nPt <= 0.05);
  fprintf('  %-44s %-10.4f\n', name, fpr);
  val_util('result', name, abs(fpr-0.05) < 0.015, sprintf('FPR %.4f', fpr));

  if scenario == 1
    t_vox = el/n_perm;
  end
end

% ---------------------------------------------------------------------
% the vectorised solve must reproduce the per-voxel pseudoinverse loop
%
% calc_GLM_voxelwise assembles Xi'*Xi and Xi'*y for every voxel at once and
% batch-solves them. That is an algebraic rearrangement, not an approximation,
% so the tolerance is the single-precision rounding floor. The reference is the
% loop the toolbox still ships, driven over every voxel. Both return the pieces
% of the statistic rather than the statistic itself, because the shrinkage of
% ResMS needs the maximum over all voxels, so the pieces are what is compared.
% ---------------------------------------------------------------------
for scenario = 1:3
  n = 30;
  C = 0.5 + 0.1*randn(n_vox, n);

  switch scenario
  case 1
    name  = 't-contrast, covariate of interest (C is permuted)';
    X     = [ones(n,1) mean(C,1)'];
    c     = [0; 1];  STAT = 'T'; eidf = 1;
    ind_X = 2; cols = 2;
  case 2
    name  = 't-contrast, covariate as nuisance (C not permuted)';
    g     = [ones(15,1); zeros(15,1)];
    X     = [g 1-g mean(C,1)'];
    c     = [1; -1; 0];  STAT = 'T'; eidf = 1;
    ind_X = 1:2; cols = 3;
  case 3
    % 3 groups, voxel-wise covariate as nuisance, F on the group effect
    name  = 'F-contrast, covariate as nuisance';
    g     = repmat((1:3)', 10, 1);
    X     = [g==1 g==2 g==3 mean(C,1)'];
    c     = [1 0; -1 1; 0 -1; 0 0];  STAT = 'F'; eidf = rank(c);
    ind_X = 1:3; cols = 4;
  end

  xX.X = double(X); xX.W = speye(n);
  xC.cols = cols;
  xCon.c = c; xCon.STAT = STAT; xCon.eidf = eidf;

  Y = single(randn(n_vox, n));

  worst = 0;
  for r = 1:10
    pr   = randperm(n);
    Pset = sparse(1:n, pr, 1, n, n);

    xXperm = xX;
    xXperm.X(:,ind_X) = Pset*xX.X(:,ind_X);

    % the reference: permute C the way the old code did, then run the loop
    Cref = C;
    if any(ismember(cols, ind_X))
      Cref = C*full(Pset);            % the n_vox x n x n multiply that was removed
    end
    Xw   = full(xXperm.W*xXperm.X);
    trRV = n - rank(xXperm.X);

    [num_r, cGc_r, ResMS_r] = ...
        calc_GLM_voxelwise_loop(Y, Xw, Cref, xC, xCon, trRV, 1, (1:n_vox)');

    % the batch gets the unpermuted C and Pset, so this also checks that the
    % column reindex reproduces the C*full(Pset) it replaced
    Cb = C;
    if any(ismember(cols, ind_X))
      [ii,jj,vv] = find(Pset);
      Cb = zeros(size(C)); Cb(:,jj) = C(:,ii).*vv(:).';
    end
    [num_b, cGc_b, ResMS_b] = ...
        calc_GLM_voxelwise_batch(Y, Xw, Cb, xC, xCon, trRV);

    rel = @(a,b) max(abs(a-b))/max(abs(b));
    worst = max([worst rel(num_b,num_r) rel(cGc_b,cGc_r) rel(ResMS_b,ResMS_r)]);
  end

  val_util('result', name, worst < 1e-5, ...
    sprintf('max rel |batch - loop| = %.2e', worst));
end

% ---------------------------------------------------------------------
% the F-statistic of a 1-df contrast must be the square of its t-statistic
%
% This pins the new voxel-wise F down against the t-path without going through
% the loop at all, so it is an independent check of the batched ESS.
% ---------------------------------------------------------------------
n = 30;
C = 0.5 + 0.1*randn(n_vox, n);
X = [ones(n,1) mean(C,1)' randn(n,1)];
xX.X = X; xX.W = speye(n);
xC.cols = 2;
Y = single(randn(n_vox, n));

xCon.c = [0;1;0]; xCon.STAT = 'T'; xCon.eidf = 1;
tt = calc_GLM_voxelwise(Y, xX, xC, xCon, ind_mask, dim, C, [], 2, 1);

xCon.STAT = 'F'; xCon.eidf = 1;
ff = calc_GLM_voxelwise(Y, xX, xC, xCon, ind_mask, dim, C, [], 2, 1);

worst = max(abs(ff(ind_mask) - tt(ind_mask).^2))/max(abs(ff(ind_mask)));
val_util('result','F of a 1-df contrast equals t^2', worst < 1e-5, ...
  sprintf('max rel |F - t^2| = %.2e', worst));

% ---------------------------------------------------------------------
% the shrinkage must use the maximum over ALL voxels, as calc_GLM does
%
% The per-voxel loop this replaced applied the shrinkage line to a scalar ResMS,
% where max(ResMS(isfinite(ResMS))) is that same scalar, so it collapsed to
% ResMS*1.001 and the voxel-wise path was in effect not shrinking at all.
% ---------------------------------------------------------------------
xCon.c = [0;1;0]; xCon.STAT = 'T'; xCon.eidf = 1;
Xw   = full(xX.W*xX.X);
trRV = n - rank(xX.X);
[num_r, cGc_r, ResMS_r] = calc_GLM_voxelwise_loop(Y, Xw, C, xC, xCon, trRV, 1, (1:n_vox)');

t_global = num_r./(eps+sqrt((ResMS_r + 1e-3*max(ResMS_r)).*cGc_r));   % what calc_GLM does
t_scalar = num_r./(eps+sqrt((ResMS_r * 1.001).*cGc_r));               % what the old code did
t_new    = calc_GLM_voxelwise(Y, xX, xC, xCon, ind_mask, dim, C, [], 2, 1);
t_new    = t_new(ind_mask);

val_util('result','shrinkage uses the maximum over all voxels', ...
  max(abs(t_new - t_global))/max(abs(t_global)) < 1e-5, ...
  sprintf('and differs from the old per-voxel one by %.1f%%', ...
          100*max(abs(t_global - t_scalar))/max(abs(t_global))));

% ---------------------------------------------------------------------
% a covariate spread over several columns must be permuted if ANY of them is
% tested, not only if the last one is
%
% With one covariate column per group, xC.cols = [1 2], a contrast on group 1
% alone tests column 1 but not column 2. The loop that used to decide this
% reassigned its flag on every pass, so the last column decided on its own and C
% was left unpermuted -- the covariate images then stayed with their original
% subjects while the design was permuted around them, which breaks the null.
% ---------------------------------------------------------------------
n  = 30;
g1 = [ones(15,1); zeros(15,1)];  g2 = 1-g1;
vc = abs(randn(n,1)) + 1;
C  = 0.5 + 0.1*randn(n_vox, n);

X = [g1.*vc g2.*vc ones(n,1)];       % one covariate column per group
xX.X = X; xX.W = speye(n);
xC.cols = [1 2];                      % the covariate lives in BOTH columns
xCon.c = [1;0;0]; xCon.STAT = 'T'; xCon.eidf = 1;   % ... but only column 1 is tested
ind_X = 1;
Y = single(randn(n_vox, n));

pr   = randperm(n);
Pset = sparse(1:n, pr, 1, n, n);
xXperm = xX;
xXperm.X(:,ind_X) = Pset*X(:,ind_X);

Xw   = full(xXperm.W*xXperm.X);
trRV = n - rank(xXperm.X);

shrink = @(num,cGc,ResMS) num./(eps+sqrt((ResMS + 1e-3*max(ResMS(isfinite(ResMS)))).*cGc));

[a,b,d] = calc_GLM_voxelwise_loop(Y, Xw, C*full(Pset), xC, xCon, trRV, 1, (1:n_vox)');
t_permuted = shrink(a,b,d);                       % what it should be
[a,b,d] = calc_GLM_voxelwise_loop(Y, Xw, C,           xC, xCon, trRV, 1, (1:n_vox)');
t_stale    = shrink(a,b,d);                       % what the old code produced

t_new = calc_GLM_voxelwise(Y, xXperm, xC, xCon, ind_mask, dim, C, Pset, ind_X, 1);
t_new = t_new(ind_mask);

val_util('result','multi-column covariate is permuted if any column is tested', ...
  max(abs(t_new - t_permuted))/max(abs(t_permuted)) < 1e-5, ...
  sprintf('the old flag left it unpermuted, which was %.0f%% off', ...
          100*max(abs(t_permuted - t_stale))/max(abs(t_permuted))));

% ---------------------------------------------------------------------
% a rank-deficient voxel must fall back to the pseudoinverse, not produce junk
%
% If the covariate is constant across subjects at some voxel, that column of Xi
% becomes collinear with the intercept and the batched solve cannot resolve it.
% Those voxels have to be handed back to pinv, which resolves them in the
% least-squares sense, exactly as before.
% ---------------------------------------------------------------------
n = 30;
C = 0.5 + 0.1*randn(n_vox, n);
C(1:50,:) = repmat(0.7, 50, n);          % 50 degenerate voxels

X = [ones(n,1) mean(C,1)'];
xX.X = X; xX.W = speye(n);
xC.cols = 2; xCon.c = [0;1]; xCon.STAT = 'T'; xCon.eidf = 1;
Y = single(randn(n_vox, n));

Xw   = full(xX.W*xX.X);
trRV = n - rank(xX.X);

[num_r, cGc_r, ResMS_r] = calc_GLM_voxelwise_loop(Y, Xw, C, xC, xCon, trRV, 1, (1:n_vox)');
[num_b, cGc_b, ResMS_b, bad] = calc_GLM_voxelwise_batch(Y, Xw, C, xC, xCon, trRV);

val_util('result','rank-deficient voxels are detected', nnz(bad) == 50, ...
  sprintf('%d of 50 degenerate voxels handed back to pinv', nnz(bad)));

% the degenerate ones must not silently pass through the batched solve: the
% pieces it produces for them have to be wrong, which is why they are handed back
val_util('result','rank-deficient voxels really do need pinv', ...
  max(abs(num_b(bad) - num_r(bad))) > 1e-3*max(abs(num_r)), ...
  'the batched solve does not resolve them');

% and away from the degenerate ones the batch must reproduce the loop
rel = @(a,b) max(abs(a(~bad)-b(~bad)))/max(abs(b(~bad)));
worst = max([rel(num_b,num_r) rel(cGc_b,cGc_r) rel(ResMS_b,ResMS_r)]);
val_util('result','the remaining voxels still match', worst < 1e-5, ...
  sprintf('max rel |batch - loop| = %.2e', worst));

% end to end: the degenerate voxels come out of calc_GLM_voxelwise via pinv
t_new = calc_GLM_voxelwise(Y, xX, xC, xCon, ind_mask, dim, C, [], 1, 1);
val_util('result','degenerate voxels survive the full path', ...
  all(isfinite(t_new(ind_mask))), 'no Inf/NaN leaks into the statistic');

% ---------------------------------------------------------------------
% cost relative to the standard (constant design) path
% ---------------------------------------------------------------------
n = 30;
X = [ones(n,1) randn(n,1)];
xX.X = X; xX.W = speye(n);
xCon.c = [0;1]; xCon.STAT = 'T';
Y = single(randn(n_vox, n));

nrep = 20;
tic; for r = 1:nrep, calc_GLM(Y, xX, xCon, ind_mask, dim); end
t_std = toc/nrep;

% the loop the batched solve replaced, for reference
C  = 0.5 + 0.1*randn(n_vox, n);
Xc = [ones(n,1) mean(C,1)'];
xCc.cols = 2; xConc.c = [0;1]; xConc.STAT = 'T';

nrep = 5;
tic; for r = 1:nrep
  calc_GLM_voxelwise_loop(Y, Xc, C, xCc, xConc, n-rank(Xc), 1, (1:n_vox)');
end
t_loop = toc/nrep;

fprintf('\n  cost per permutation (%d elements, n = %d):\n', n_vox, n);
fprintf('    standard design (constant)    : %8.4f s\n', t_std);
fprintf('    voxel-wise, per-voxel pinv    : %8.4f s   (%.0fx slower)\n', ...
        t_loop, t_loop/t_std);
fprintf('    voxel-wise, batched solve     : %8.4f s   (%.0fx slower)\n', ...
        t_vox, t_vox/t_std);
fprintf('    -> extrapolated to 400k elements, 5000 permutations:\n');
fprintf('       standard          : %6.1f h\n', t_std*(400000/n_vox)*5000/3600);
fprintf('       per-voxel pinv    : %6.1f h\n', t_loop*(400000/n_vox)*5000/3600);
fprintf('       batched solve     : %6.1f h\n', t_vox*(400000/n_vox)*5000/3600);

val_util('result','the batched solve is faster than the pinv loop', ...
  t_loop/t_vox > 2, sprintf('%.1fx faster than per-voxel pinv', t_loop/t_vox));

val_util('summary');
