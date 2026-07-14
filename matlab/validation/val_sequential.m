function val_sequential
% Does sequential stopping reach the same answer as running every permutation?
%
% Instead of always running n_perm permutations, the loop stops once the observed
% GLOBAL maximum has been exceeded n_exceed_stop times by the permutation maxima
% (Besag & Clifford, 1991; the negative binomial method of Winkler et al., 2016).
%
% The reason this is the right thing to watch is that no element rests on fewer
% exceedances than the global maximum does: an element with a smaller statistic is
% exceeded at least as often by the permutation maxima. So once the global maximum
% has been exceeded n_exceed_stop times, every corrected p-value in the image is
% known to a relative standard error of 1/sqrt(n_exceed_stop) or better.
%
% The saving is entirely one-sided, and that is the point. An image with nothing
% in it has a largest value that is an ordinary draw from the very distribution it
% is compared against, so it is exceeded about every other permutation and the
% target is reached at once. An image with a real effect has a largest value that
% is almost never exceeded, never reaches the target, and runs the full n_perm.
%
% What has to be established is that the early stop does not change the answer:
%   * the FWE decision must agree with the decision the full run would have made
%   * the false-positive rate must still be alpha
%   * an image with a real effect must NOT be stopped early
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Sequential stopping');

rng(11);

n_perm        = 5000;    % what a full run would do
n_perm_min    = 500;     % the floor, as in tfce_estimate_stat
n_exceed_stop = 20;      % as in tfce_estimate_stat
n_sigma_stop  = 3;       % as in tfce_estimate_stat
alpha         = 0.05;
n_sim         = 400;     % simulated analyses
n_elem        = 200;     % independent elements per analysis

% the stopping rule of tfce_estimate_stat, applied to a stream of permutation
% maxima: enough exceedances, AND the estimate is n_sigma_stop standard errors
% clear of alpha
  function m = stop_at(perm_max, obs_max)
    e     = cumsum(perm_max >= obs_max);
    mm    = 1:numel(perm_max);
    p_hat = e./mm;
    se    = sqrt(p_hat.*(1-p_hat)./mm);
    ok    = (mm >= n_perm_min) & (e >= n_exceed_stop) & ...
            ((p_hat - n_sigma_stop*se) > alpha);
    m     = find(ok, 1);
    if isempty(m), m = numel(perm_max); end
  end

% ---------------------------------------------------------------------
% Each simulated analysis needs a permutation distribution of the maximum and an
% observed maximum. The maximum of many correlated positive statistics is well
% described by a Gumbel-like law; what matters here is only that the observed
% maximum is a draw from the SAME law as the permuted ones under the null, and
% from a shifted one when there is an effect.
% ---------------------------------------------------------------------
% -log(-log(u)) with u uniform is a standard Gumbel draw, so this needs nothing
% beyond base MATLAB
draw_max = @(n, shift) max(-log(-log(rand(n_elem, n))), [], 1) + shift;

for scenario = 1:2

  if scenario == 1
    shift = 0;                       % null: nothing in the image
    name  = 'null image';
  else
    shift = 8;                       % a strong effect
    name  = 'image with a real effect';
  end

  stopped_at = zeros(n_sim,1);
  agree      = true(n_sim,1);
  rej_seq    = false(n_sim,1);
  rej_full   = false(n_sim,1);

  for s = 1:n_sim
    perm_max = draw_max(n_perm, 0);        % the null distribution of the maximum
    obs_max  = draw_max(1, shift);         % the observed maximum

    hit = stop_at(perm_max, obs_max);      % where the rule would have stopped
    stopped_at(s) = hit;

    % the FWE decision from the truncated run, and from the full one
    rej_seq(s)  = (sum(perm_max(1:hit) >= obs_max)/hit)     <= alpha;
    rej_full(s) = (sum(perm_max >= obs_max)/n_perm)         <= alpha;
    agree(s)    = rej_seq(s) == rej_full(s);
  end

  fprintf('\n  %s\n', name);
  fprintf('    permutations run : median %d, worst case %d (of %d)\n', ...
          round(median(stopped_at)), max(stopped_at), n_perm);
  fprintf('    saving           : %.1fx\n', n_perm/mean(stopped_at));
  fprintf('    FWE rejections   : sequential %.3f, full run %.3f\n', ...
          mean(rej_seq), mean(rej_full));

  if scenario == 1
    % a null image must be stopped as early as the floor allows
    val_util('result','null image stops at the floor', ...
      median(stopped_at) <= n_perm_min + 10, ...
      sprintf('median %d permutations instead of %d (%.0fx)', ...
              round(median(stopped_at)), n_perm, n_perm/mean(stopped_at)));

    % and the false-positive rate must be unchanged
    val_util('result','false-positive rate is still alpha', ...
      abs(mean(rej_seq) - alpha) < 0.02, ...
      sprintf('%.3f sequential vs %.3f for the full run (nominal %.2f)', ...
              mean(rej_seq), mean(rej_full), alpha));
  else
    % a real effect must never be stopped early: its maximum is not exceeded
    val_util('result','a real effect is never stopped early', ...
      all(stopped_at == n_perm), ...
      sprintf('all %d analyses ran the full %d permutations', n_sim, n_perm));

    val_util('result','a real effect is still detected', ...
      mean(rej_seq) > 0.99, ...
      sprintf('rejected in %.1f%% of analyses', 100*mean(rej_seq)));
  end

  % the decision must be the one the full run would have made
  val_util('result', sprintf('%s: decision matches the full run', name), ...
    mean(agree) > 0.98, ...
    sprintf('%.1f%% agreement over %d analyses', 100*mean(agree), n_sim));
end

% ---------------------------------------------------------------------
% the borderline case is the one that can actually go wrong
%
% An image whose true corrected p-value sits right at alpha is where a truncated
% run and a full run would be most likely to disagree, because both estimates are
% noisy and both straddle the threshold. Counting exceedances alone does not
% protect against this at all: the target number of exceedances is reached quickly
% even at p = alpha, and the decision would then rest on an estimate far too
% coarse to place it on one side or the other. Requiring the estimate to be
% n_sigma_stop standard errors CLEAR of alpha is what protects against it -- such
% an image can never satisfy that, so it is never cut short.
% ---------------------------------------------------------------------
n_sim  = 2000;
agree  = true(n_sim,1);
early  = false(n_sim,1);
seq_p  = zeros(n_sim,1);
full_p = zeros(n_sim,1);

% an observed maximum whose true corrected p-value sits exactly at alpha
ref     = sort(draw_max(200000, 0), 'descend');
obs_max = ref(round(alpha*numel(ref)));

for s = 1:n_sim
  perm_max = draw_max(n_perm, 0);

  hit      = stop_at(perm_max, obs_max);
  early(s) = hit < n_perm;

  seq_p(s)  = sum(perm_max(1:hit) >= obs_max)/hit;
  full_p(s) = sum(perm_max >= obs_max)/n_perm;
  agree(s)  = (seq_p(s) <= alpha) == (full_p(s) <= alpha);
end

fprintf('\n  images sitting exactly at alpha = %.2f\n', alpha);
fprintf('    stopped early    : %.1f%% of %d analyses\n', 100*mean(early), n_sim);
fprintf('    decision agrees  : %.1f%%\n', 100*mean(agree));

% The resampling risk: how often does an image sitting exactly at alpha get cut
% short at all? A n_sigma_stop = 3 margin is a one-sided 0.1% chance on a single
% look, and the rule looks after every permutation, so a little more than that
% leaks through. It stays small, and the check below shows that even those runs
% reach the same decision.
val_util('result','borderline images are almost never cut short', ...
  mean(early) < 0.02, ...
  sprintf('%.1f%% of %d stopped early (a %d-SE margin, tested at every permutation)', ...
          100*mean(early), n_sim, n_sigma_stop));

val_util('result','borderline images: decision matches the full run', ...
  mean(agree) > 0.99, ...
  sprintf('%.1f%% agreement at exactly the alpha boundary', 100*mean(agree)));

val_util('summary');
end
