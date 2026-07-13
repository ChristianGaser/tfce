function val_nuisance
% Calibration of the three strategies for nuisance variables.
%
% Draper-Stoneman permutes the effects of interest in the design, Freedman-Lane
% permutes the data after residualising it against the nuisance variables, and
% Smith orthogonalises the effects of interest against the nuisance space before
% permuting them (Winkler et al., 2014). All three are exact only under
% conditions that differ, and Draper-Stoneman is known to behave poorly when the
% nuisance variables are correlated with the effect of interest.
%
% The permutations are built exactly as tfce_estimate_stat builds them, so this
% is a check of the implementation, not of the published methods.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Nuisance variables: Draper-Stoneman / Freedman-Lane / Smith');

rng(2);
n1 = 15; n2 = 15; n = n1+n2;
n_vox = 3000; n_perm = 500;
g = [ones(n1,1); zeros(n2,1)];

scenarios = {
  'no nuisance',                     0.0
  'nuisance uncorrelated with group', 0.0
  'nuisance correlated with group',   2.0
};

fprintf('\n  uncorrected false-positive rate at alpha = 0.05 (nominal 0.050)\n');
fprintf('  %-34s %-12s %-14s %-10s\n','scenario','Draper-Ston.','Freedman-Lane','Smith');
fprintf('  %s\n', repmat('-',1,74));

for s = 1:size(scenarios,1)
  name = scenarios{s,1};
  rho  = scenarios{s,2};

  if s == 1
    X = [g 1-g];  c = [1;-1];  ind_X = 1:2;  ind_Z = [];
  else
    cov1 = rho*g + randn(n,1);              % correlated with group if rho > 0
    X = [g 1-g cov1];  c = [1;-1;0];  ind_X = 1:2;  ind_Z = 3;
  end

  % null data, but the nuisance covariate has a STRONG real effect
  Y = randn(n, n_vox);
  if ~isempty(ind_Z)
    Y = Y + 3*X(:,ind_Z)*ones(1,n_vox);
  end

  % Guttman partitioning, exactly as get_nuisance_method does
  if isempty(ind_Z)
    Rz = eye(n);
  else
    Z  = X(:,ind_Z);
    Rz = eye(n) - Z*pinv(Z);
  end

  t0  = glm_t(Y, X, c);
  fpr = zeros(1,3);

  for method = 0:2
    cnt = zeros(1,n_vox);
    for p = 1:n_perm
      pr   = randperm(n);
      Pset = sparse(1:n, pr, 1, n, n);       % permutation matrix

      switch method
        case 0   % Draper-Stoneman: permute the effects of interest in X
          Xp = X; Xp(:,ind_X) = Pset*X(:,ind_X);
          tp = glm_t(Y, Xp, c);
        case 1   % Freedman-Lane: permute the residualised data, X unchanged
          Yp = (Pset'*Rz)' * Y;              % = Rz' * Pset * Y  (rows = data)
          tp = glm_t(Yp, X, c);
        case 2   % Smith: orthogonalise, then permute
          Xp = X; Xp(:,ind_X) = Pset*Rz*X(:,ind_X);
          tp = glm_t(Y, Xp, c);
      end
      cnt = cnt + (tp >= t0);
    end
    nPt = cnt/n_perm;
    fpr(method+1) = mean(nPt <= 0.05);
  end

  fprintf('  %-34s %-12.4f %-14.4f %-10.4f\n', name, fpr(1), fpr(2), fpr(3));

  tol = 0.015;
  val_util('result', [name ' - Draper-Stoneman'],  abs(fpr(1)-0.05) < tol, sprintf('FPR %.4f', fpr(1)));
  val_util('result', [name ' - Freedman-Lane'],    abs(fpr(2)-0.05) < tol, sprintf('FPR %.4f', fpr(2)));
  val_util('result', [name ' - Smith'],            abs(fpr(3)-0.05) < tol, sprintf('FPR %.4f', fpr(3)));
end

val_util('summary');

%---------------------------------------------------------------
function t = glm_t(Y, X, c)
pKX = pinv(X); n = size(X,1);
B = pKX*Y;
R = Y - X*B;
ResMS = sum(R.^2,1)/(n - rank(X));
t = (c'*B)./sqrt(ResMS*(c'*(pKX*pKX')*c));
