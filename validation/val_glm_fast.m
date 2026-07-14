function val_glm_fast
% Is the accelerated GLM identical to the unaccelerated one, and how much faster?
%
% calc_GLM_fast obtains the permuted statistic from the n_vox x rank(X) matrix
% Beta alone. It never forms the permuted data matrix and never forms the
% residual matrix, because
%
%   ResSS = ||R*y||^2 = ||y||^2 - Beta'*(X'*X)*Beta
%
% and, for Freedman-Lane, because the nuisance columns Z are part of the full
% design X, so the full-model residual forming matrix already annihilates them
% (R*Rz = R) and Rz drops out of the residuals:
%
%   ResSS = ||R*Rz*Pset*y||^2 = ||R*Pset*y||^2 = ||Pset*y||^2 - Beta'*(X'*X)*Beta
%
% Since Pset only permutes or sign-flips the rows it touches, ||Pset*y||^2 is the
% same for every permutation. That removes the n_vox x n_data x n_data multiply
% Y*(Pset'*Rz) from the permutation loop entirely.
%
% This is an algebraic identity, not an approximation, so the two paths must
% agree to single-precision rounding. The checks below establish that on t- and
% F-contrasts for all three nuisance methods.
%
% R*Rz = R only holds if Z lies in the space spanned by the design that is
% actually fitted. calc_GLM fits W*X on data whitened at load time, while Rz is
% built from the unwhitened Z, so the identity breaks for a non-identity
% whitening matrix. The last check establishes that prepare_GLM_fast detects
% this and refuses the fast path, so the loop falls back to calc_GLM.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Accelerated GLM');
val_util('extract','calc_GLM','calc_GLM_fast','prepare_GLM_fast','make_Pset', ...
         'permuted_design');

rng(1);
tol = 1e-5;   % Beta is single, so this is the rounding floor, not slack

methods = {'Draper-Stoneman','Freedman-Lane','Smith'};

% ---------------------------------------------------------------------
% the two paths must produce the same statistic
% ---------------------------------------------------------------------
for STAT = 'TF'
  for nuisance_method = 0:2
    [Y, xX, xCon, Rz, ind_X, n_cond, ind_label] = make_design(60, 3000, STAT);
    n = size(xX.X,1);

    pre = prepare_GLM_fast(Y, xX, Rz, nuisance_method, ...
        make_Pset(randperm(n), n_cond, n, n, ind_label), ...
        make_Pset(randperm(n), n_cond, n, n, ind_label));

    worst = 0;
    accepted = ~isempty(pre);
    if accepted
      for r = 1:20
        Pset  = make_Pset(randperm(n), n_cond, n, n, ind_label);
        xXperm = permute_design(xX, Pset, Rz, ind_X, nuisance_method);

        if nuisance_method == 1
          t_ref = calc_GLM(Y*(Pset'*Rz), xXperm, xCon, (1:size(Y,1))', [size(Y,1) 1]);
        else
          t_ref = calc_GLM(Y, xXperm, xCon, (1:size(Y,1))', [size(Y,1) 1]);
        end
        t_new = calc_GLM_fast(Y, pre, xXperm, xCon, (1:size(Y,1))', [size(Y,1) 1], ...
            Pset, nuisance_method);

        worst = max(worst, max(abs(t_new(:)-t_ref(:)))/max(abs(t_ref(:))));
      end
    end

    val_util('result', ...
      sprintf('%s-contrast, %s', STAT, methods{nuisance_method+1}), ...
      accepted && worst < tol, ...
      sprintf('max rel |t_fast - t_ref| = %.2e', worst));
  end
end

% ---------------------------------------------------------------------
% the fast path must refuse itself when R*Rz = R does not hold
% ---------------------------------------------------------------------
[Y, xX, ~, Rz, ~, n_cond, ind_label] = make_design(60, 500, 'T');
n  = size(xX.X,1);
P0 = make_Pset(randperm(n), n_cond, n, n, ind_label);
P1 = make_Pset(randperm(n), n_cond, n, n, ind_label);

val_util('result','identity whitening: fast path is taken', ...
  ~isempty(prepare_GLM_fast(Y, xX, Rz, 1, P0, P1)), 'W = I');

w = ones(n,1); w(1:2:end) = 2;
xXw = xX; xXw.W = spdiags(w, 0, n, n);

val_util('result','non-identity whitening: Freedman-Lane refuses', ...
  isempty(prepare_GLM_fast(Y, xXw, Rz, 1, P0, P1)), ...
  'R*Rz ~= R, so the loop falls back to calc_GLM');

val_util('result','non-identity whitening: Draper-Stoneman still fast', ...
  ~isempty(prepare_GLM_fast(Y, xXw, Rz, 0, P0, P1)), ...
  'does not rely on R*Rz = R');

% ---------------------------------------------------------------------
% what it buys
% ---------------------------------------------------------------------
fprintf('\n  cost per permutation (400k elements):\n');
fprintf('    %-28s %-11s %-11s %s\n','design','calc_GLM','fast','speedup');
fprintf('  %s\n', repmat('-',1,62));

speedup = zeros(1,2);
cases   = {{'n = 100, 2 nuisance', 100}, {'n = 250, 4 nuisance', 250}};

for k = 1:numel(cases)
  n = cases{k}{2};
  [Y, xX, xCon, Rz, ind_X, n_cond, ind_label] = make_design(n, 400000, 'T');
  ind_mask = (1:size(Y,1))'; dim = [size(Y,1) 1];

  Psets = cell(1,10);
  for r = 1:10, Psets{r} = make_Pset(randperm(n), n_cond, n, n, ind_label); end
  pre = prepare_GLM_fast(Y, xX, Rz, 1, Psets{1}, Psets{2});

  tic
  for r = 1:10
    xXperm = permute_design(xX, Psets{r}, Rz, ind_X, 1);
    calc_GLM(Y*(Psets{r}'*Rz), xXperm, xCon, ind_mask, dim);
  end
  t_ref = toc/10;

  tic
  for r = 1:10
    xXperm = permute_design(xX, Psets{r}, Rz, ind_X, 1);
    calc_GLM_fast(Y, pre, xXperm, xCon, ind_mask, dim, Psets{r}, 1);
  end
  t_new = toc/10;

  speedup(k) = t_ref/t_new;
  fprintf('    %-28s %-11s %-11s %.1fx\n', cases{k}{1}, ...
    sprintf('%.3f s', t_ref), sprintf('%.3f s', t_new), speedup(k));
end

fprintf('    -> GLM time for 5000 permutations, Freedman-Lane, n = 250:\n');
fprintf('       calc_GLM  : %5.1f min\n', 5000*t_ref/60);
fprintf('       fast      : %5.1f min\n', 5000*t_new/60);

val_util('result','the fast path is actually faster', all(speedup > 1.5), ...
  sprintf('%.1fx / %.1fx on the GLM', speedup(1), speedup(2)));

val_util('summary');


%---------------------------------------------------------------
function xXperm = permute_design(xX, Pset, Rz, ind_X, nuisance_method)
% apply a permutation to the design, through the very function the loop uses
%
% n_exch_blocks = 1 keeps permuted_design out of its interaction-design branch,
% which these designs do not have.

xXperm   = xX;
xXperm.X = permuted_design(xX.X, Pset, Rz, ind_X, nuisance_method, 1, 0, true);


%---------------------------------------------------------------
function [Y, xX, xCon, Rz, ind_X, n_cond, ind_label] = make_design(n, n_vox, STAT)
% a second-level design with nuisance variables, in the layout the toolbox uses

if STAT == 'F'
  n_grp = 3;
  g  = repmat((1:n_grp)', ceil(n/n_grp), 1); g = g(1:n);
  Xi = zeros(n, n_grp);
  for k = 1:n_grp, Xi(g==k,k) = 1; end
else
  Xi = zeros(n,2);
  Xi(1:round(n/2),1)     = 1;
  Xi(round(n/2)+1:end,2) = 1;
end

Z = randn(n, 2);                 % age, TIV
if n >= 250, Z = [Z randn(n,2)]; end

X     = [Xi Z];
p     = size(Xi,2);
r     = size(X,2);
ind_X = 1:p;

xX.X  = X;
xX.W  = speye(n);
xX.iH = 1:p;
xX.iC = p+1:r;
xX.iB = [];
xX.iG = [];

if STAT == 'T'
  xCon.c    = [1; -1; zeros(r-2,1)];
  xCon.STAT = 'T';
  xCon.eidf = 1;
else
  C = zeros(r, p-1);
  for k = 1:p-1, C(k,k) = 1; C(k+1,k) = -1; end
  xCon.c    = C;
  xCon.STAT = 'F';
  xCon.eidf = rank(C);
end

Rz = eye(n) - Z*pinv(Z);

Y = single(randn(n_vox, n));
Y = Y + single(0.3 * X(:,1)');   % a little signal, so ResSS is not pure noise

n_cond    = p;                   % > 1, i.e. permutation rather than sign-flipping
ind_label = (1:n)';
