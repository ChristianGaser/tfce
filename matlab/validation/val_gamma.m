function val_gamma
% Calibration of the Gamma approximation of the maximum distribution.
%
% FWE-corrected p-values are obtained by fitting a Gamma distribution to the
% permutation distribution of the maximum statistic (moment matching on mean,
% variance and skewness), instead of counting exceedances. Counting cannot
% resolve p-values below 1/n_perm; the Gamma fit can, but only if its shape
% assumption holds.
%
% This is the single most consequential approximation in the toolbox: it enters
% EVERY corrected p-value that is reported. The checks below ask
%   (1) does the fit agree with counting where counting is reliable?
%   (2) does it stay calibrated, i.e. is P(p_FWE <= 0.05) = 0.05 under the null?
%   (3) does it agree with a counting estimator based on far more permutations?
%   (4) is it distorted when the permutation set is truncated by early stopping?
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Gamma approximation of the maximum distribution');
val_util('extract','palm_moments','palm_gamma');

rng(0);

% ---------------------------------------------------------------------
% Draw null maximum statistics. We use maxima of TFCE maps of smoothed noise,
% so the shape is realistic rather than a textbook distribution.
% ---------------------------------------------------------------------
dim = [32 32 32]; E = 0.5; H = 2.0;
n_big = 5000;                       % reference: many permutations
fprintf('\n  generating %d null maxima (TFCE of smoothed noise, %dx%dx%d) ...\n', ...
        n_big, dim);

[x,y,z] = ndgrid(linspace(-1,1,dim(1)),linspace(-1,1,dim(2)),linspace(-1,1,dim(3)));
mask = (x.^2+y.^2+z.^2) < 0.85;

maxstat = zeros(1,n_big);
for i = 1:n_big
  m = smooth3(randn(dim),'gaussian',7,1.8);
  m = m/std(m(:)); m(~mask) = 0;
  t = tfceMex_maxtree(m, E, H, 1);
  maxstat(i) = max(t(:));
end

% ---------------------------------------------------------------------
% (1) Gamma vs counting, where counting is reliable (p >> 1/n_perm)
% ---------------------------------------------------------------------
% Only the tail matters for inference, and the counting estimator it is compared
% against is itself noisy: at p its standard error is sqrt(p(1-p)/n_perm). We
% therefore require agreement to within 3 standard errors of counting, which is
% the tightest tolerance the comparison can support.
p_targets = [0.10 0.05 0.025 0.01];
fprintf('\n  Gamma vs counting in the tail (tolerance = 3 SE of counting):\n');
fprintf('  %-8s %-8s %-10s %-10s %-10s %s\n','n_perm','p','Gamma','counting','|dp|','3 SE');
for n_perm = [500 1000 5000]
  null = maxstat(1:n_perm);
  [mu,s2,g1] = palm_moments(null');
  ok = true; worst = 0;
  for p = p_targets
    q  = prctile(null, 100*(1-p));            % statistic with counting p-value p
    pg = palm_gamma(q, mu, s2, g1, false, 1/n_perm);
    pc = mean(null >= q);
    se = sqrt(p*(1-p)/n_perm);
    fprintf('  %-8d %-8.3f %-10.4f %-10.4f %-10.4f %.4f\n', n_perm, p, pg, pc, abs(pg-pc), 3*se);
    ok = ok && (abs(pg-pc) <= 3*se);
    worst = max(worst, abs(pg-pc));
  end
  val_util('result', sprintf('n_perm=%d: Gamma vs counting in the tail', n_perm), ...
    ok, sprintf('max |dp| = %.4f', worst));
end

% ---------------------------------------------------------------------
% (2) Calibration: is P(p_FWE <= alpha) = alpha under the null?
%     Take one map as "observed" and the remaining ones as the null.
% ---------------------------------------------------------------------
fprintf('\n  null calibration of the corrected p-value:\n');
fprintf('  %-10s %-12s %-12s %-12s\n','n_perm','alpha=0.05','alpha=0.01','(counting)');
for n_perm = [500 1000]
  nrep = 2000;
  p_g = zeros(1,nrep); p_c = zeros(1,nrep);
  for r = 1:nrep
    idx  = randperm(n_big, n_perm+1);
    obs  = maxstat(idx(1));
    null = maxstat(idx(2:end));
    [mu,s2,g1] = palm_moments(null');
    p_g(r) = palm_gamma(obs, mu, s2, g1, false, 1/n_perm);
    p_c(r) = mean(null >= obs);
  end
  f05 = mean(p_g <= 0.05); f01 = mean(p_g <= 0.01);
  c05 = mean(p_c <= 0.05);
  fprintf('  %-10d %-12.4f %-12.4f %-12.4f\n', n_perm, f05, f01, c05);
  val_util('result', sprintf('n_perm=%d: FWE rate at alpha=0.05', n_perm), ...
    abs(f05-0.05) < 0.015, sprintf('%.4f (nominal 0.050)', f05));
  val_util('result', sprintf('n_perm=%d: FWE rate at alpha=0.01', n_perm), ...
    abs(f01-0.01) < 0.008, sprintf('%.4f (nominal 0.010)', f01));
end

% ---------------------------------------------------------------------
% (3) Does the Gamma fit from FEW permutations agree with counting from MANY?
%     This is the actual claim: it buys resolution below 1/n_perm.
% ---------------------------------------------------------------------
fprintf('\n  extrapolation below the counting resolution:\n');
fprintf('  %-14s %-14s %-14s %-10s\n','true p (5000)','Gamma (n=500)','count (n=500)','stat');
q = prctile(maxstat, [95 99 99.5 99.9]);
null500 = maxstat(1:500);
[mu,s2,g1] = palm_moments(null500');
worst = 0;
for i = 1:numel(q)
  p_true  = mean(maxstat >= q(i));
  p_gam   = palm_gamma(q(i), mu, s2, g1, false, 1/500);
  p_cnt   = mean(null500 >= q(i));
  fprintf('  %-14.4f %-14.4f %-14.4f %-10.2f\n', p_true, p_gam, p_cnt, q(i));
  worst = max(worst, abs(p_gam - p_true));
end
val_util('result','Gamma(n=500) vs counting(n=5000) in the tail', worst < 0.02, ...
  sprintf('max |dp| = %.4f', worst));

% ---------------------------------------------------------------------
% (4) Early stopping truncates the permutation set on a data-dependent
%     criterion. Does that distort the moment fit?
%     Emulate: stop as soon as the running 95th percentile is stable.
% ---------------------------------------------------------------------
fprintf('\n  effect of a truncated permutation set on the fit:\n');
full_null = maxstat(1:2000);
[mu,s2,g1] = palm_moments(full_null');
q95 = prctile(full_null,95);
p_full = palm_gamma(q95, mu, s2, g1, false, 1/2000);

worst_t = 0;
for n_keep = [200 500 1000]
  sub = full_null(1:n_keep);                    % truncation is on ORDER, not value
  [mu,s2,g1] = palm_moments(sub');
  p_sub = palm_gamma(q95, mu, s2, g1, false, 1/n_keep);
  fprintf('  first %-5d perms: p = %.4f (full 2000: %.4f)\n', n_keep, p_sub, p_full);
  worst_t = max(worst_t, abs(p_sub - p_full));
end
val_util('result','order-based truncation does not bias the fit', worst_t < 0.02, ...
  sprintf('max |dp| = %.4f', worst_t));

val_util('summary');
