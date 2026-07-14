function val_pareto
% Does the Generalised Pareto tail approximation resolve P-values below 1/n_perm?
%
% Counting exceedances cannot report a P-value smaller than 1/n_perm, and that
% floor -- not the FWE-corrected null, which the Gamma fit already handles -- is
% what forces a permutation test to run many thousands of permutations. The
% toolbox therefore fits a Generalised Pareto distribution to the tail of the
% permutation distribution of every element and reads the uncorrected P-value off
% that fit. Those P-values also feed the FDR correction, so the floor caps that
% too.
%
% The honest question is not whether the fit is exact -- extrapolating a tail from
% a hundred points cannot be -- but whether it beats what it replaces. So the
% checks below build a real permutation distribution from many permutations, show
% the fit only the first few, and compare BOTH the fit and plain counting against
% the truth. The fit has to win, or at least not lose, everywhere.
%
% Two properties matter as much as accuracy:
%
%   * It must never return zero. Counting does, constantly, as soon as the true P
%     drops below 1/n_perm, and a zero P-value becomes Inf in the -log10 map that
%     the toolbox writes out. Probability-weighted moments put the shape k above
%     zero for about half the elements by chance, which gives the fitted
%     distribution a finite upper end point; where the observed statistic lies
%     beyond it the fit would return zero as well. Such a fit is contradicted by
%     the observation itself, so it is rejected and the exponential limit, which
%     has infinite support, answers instead.
%
%   * It must leave alone the elements that counting already resolves.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Pareto tail approximation of uncorrected P-values');
val_util('extract','tail_init','tail_update','tail_pvalue', ...
         'gpd_pwm','gpd_sf','gpd_anderson');

rng(7);

n         = 40;      % subjects
n_vox     = 500;     % independent elements
n_ref     = 50000;   % permutations behind the reference
n_exc_min = 25;      % as in tfce_estimate_stat
K         = 100;     % as in tfce_estimate_stat (n_tail)

% ---------------------------------------------------------------------
% a genuine permutation distribution: one-sample t, by sign-flipping
% ---------------------------------------------------------------------
Y   = randn(n_vox, n);
REF = zeros(n_ref, n_vox, 'single');

for p = 1:n_ref
  s  = sign(randn(1,n));
  Ys = Y .* s;
  REF(p,:) = single((mean(Ys,2)./(std(Ys,0,2)/sqrt(n)))');
end

REFs = sort(REF, 1, 'descend');
qs   = [3e-3 1e-3 3e-4 1e-4];

% ---------------------------------------------------------------------
% the buffer must really hold the largest values of the stream
% ---------------------------------------------------------------------
tb = tail_init(K, n_vox);
for p = 1:1000, tb = tail_update(tb, REF(p,:)); end

want = sort(REF(1:1000,:), 1, 'descend');
val_util('result','the tail buffer keeps the true top-K', ...
  isequal(want(1:K,:), sort(tb.buf, 1, 'descend')), ...
  sprintf('K = %d over 1000 permutations, all %d elements', K, n_vox));

% ---------------------------------------------------------------------
% accuracy, against a reference count over n_ref permutations
% ---------------------------------------------------------------------
for n_acc = [1000 5000]

  tb = tail_init(K, n_vox);
  for p = 1:n_acc, tb = tail_update(tb, REF(p,:)); end

  fprintf('\n  P-hat / P-true from %d permutations (reference: a count over %d)\n', ...
          n_acc, n_ref);
  fprintf('    %-9s %-10s %7s %7s %7s %7s %7s %8s %9s\n', ...
          'true P','estimator','5%','25%','50%','75%','95%','P = 0','within 2x');
  fprintf('  %s\n', repmat('-',1,82));

  for q = qs
    x   = double(REFs(round(q*n_ref), :));
    cnt = sum(double(REF(1:n_acc,:)) >= x, 1);

    p_par = tail_pvalue(tb, single(x), cnt, n_acc, n_exc_min);
    p_cnt = cnt/n_acc;

    r_par = p_par/q;
    r_cnt = p_cnt/q;

    % how often does the estimate land within a factor of two of the truth? This
    % is the measure to compare the two on, rather than the spread: an estimator
    % that returns zero every time has no spread at all and is still useless, and
    % counting does exactly that once the true P drops below 1/n_perm.
    hit = @(r) 100*mean(r >= 0.5 & r <= 2);

    fprintf('    %-9.0e %-10s %7.2f %7.2f %7.2f %7.2f %7.2f %8d %8.0f%%\n', ...
            q, 'Pareto',   prctile(r_par,[5 25 50 75 95]), nnz(p_par==0), hit(r_par));
    fprintf('    %-9s %-10s %7.2f %7.2f %7.2f %7.2f %7.2f %8d %8.0f%%\n', ...
            '', 'counting', prctile(r_cnt,[5 25 50 75 95]), nnz(p_cnt==0), hit(r_cnt));

    % the fit must be unbiased in the middle. The reference is itself a count, so
    % allow for its own noise: 3 SE of q*n_ref exceedances, plus 25% for the fit.
    se   = sqrt(q*(1-q)/n_ref)/q;
    band = 1.25 + 3*se;
    med  = median(r_par);

    val_util('result', sprintf('%5d perms, P = %.0e: fit is unbiased', n_acc, q), ...
      med > 1/band && med < band, ...
      sprintf('median P-hat/P-true = %.2f (band %.2f)', med, band));

    % and it must not be less accurate than the counting it overrides
    val_util('result', sprintf('%5d perms, P = %.0e: no worse than counting', n_acc, q), ...
      hit(r_par) >= hit(r_cnt) - 5, ...
      sprintf('within a factor of 2: %.0f%% vs %.0f%% counting', ...
              hit(r_par), hit(r_cnt)));

    % a floor on the accuracy, so that a regression in the fit is caught rather
    % than merely being no worse than a counting estimate that is itself hopeless
    floors = containers.Map({'1000_3e-03','1000_1e-03','1000_3e-04','1000_1e-04', ...
                             '5000_3e-03','5000_1e-03','5000_3e-04','5000_1e-04'}, ...
                            {82, 67, 50, 36, 93, 87, 73, 54});
    key    = sprintf('%d_%.0e', n_acc, q);
    floor_ = floors(key);

    val_util('result', sprintf('%5d perms, P = %.0e: meets the accuracy floor', n_acc, q), ...
      hit(r_par) >= floor_, ...
      sprintf('%.0f%% within a factor of 2 (floor %d%%)', hit(r_par), floor_));

    % never zero, whatever counting does
    val_util('result', sprintf('%5d perms, P = %.0e: never returns zero', n_acc, q), ...
      ~any(p_par == 0) && all(isfinite(p_par)), ...
      sprintf('counting returned zero for %d of %d elements', nnz(p_cnt==0), n_vox));
  end
end

% ---------------------------------------------------------------------
% the premise of pooling: the elements share a shape, and differ in scale
%
% If that is true, then the shape estimated per element is a noisy view of one
% common number: the individual estimates must scatter widely, while their median
% must be reproducible from one half of the image to the other. If instead the
% elements genuinely had different shapes, the two halves would disagree and
% pooling would be the wrong thing to do.
% ---------------------------------------------------------------------
tb = tail_init(K, n_vox);
for p = 1:1000, tb = tail_update(tb, REF(p,:)); end

Gall = double(sort(tb.buf, 1, 'descend'));
m    = 60;
u    = (Gall(m,:) + Gall(m+1,:))/2;
kv   = gpd_pwm(flipud(Gall(1:m,:)) - u);
kv   = kv(isfinite(kv));

h1 = median(kv(1:2:end));
h2 = median(kv(2:2:end));

val_util('result','the shape is common to the elements, the noise is not', ...
  abs(h1-h2) < 0.1 && std(kv) > 3*abs(h1-h2), ...
  sprintf('per-element k scatters by %.2f, but its median agrees to %.3f across halves', ...
          std(kv), abs(h1-h2)));

% ---------------------------------------------------------------------
% elements that counting already resolves must be left exactly alone
% ---------------------------------------------------------------------

x     = double(REFs(round(0.2*n_ref), :));      % reference P = 0.2, ~200 exceedances
cnt   = sum(double(REF(1:1000,:)) >= x, 1);
p_par = tail_pvalue(tb, single(x), cnt, 1000, n_exc_min);

val_util('result','well-counted elements are not touched by the fit', ...
  isequal(p_par, cnt/1000), ...
  sprintf('all %d had >= %d exceedances and kept their count', n_vox, n_exc_min));

% ---------------------------------------------------------------------
% the fit must join the count continuously where it takes over
%
% At the threshold u the fitted survival term is 1 by construction, so the fitted
% P-value there is exactly m/n_draws. A discontinuity would show up as a jump
% across the point where an element stops being counted and starts being fitted.
% ---------------------------------------------------------------------
worst = 0;
for q = [0.024 0.026]                           % straddles n_exc_min/1000 = 0.025
  x   = double(REFs(round(q*n_ref), :));
  cnt = sum(double(REF(1:1000,:)) >= x, 1);
  pp  = tail_pvalue(tb, single(x), cnt, 1000, n_exc_min);
  worst = max(worst, median(abs(pp - cnt/1000))/q);
end

val_util('result','the fit joins the count where it takes over', worst < 0.25, ...
  sprintf('median jump across the switch = %.1f%% of P', 100*worst));

val_util('summary');
