function val_calibration
% Can the calibration check of the permutation null actually answer its question?
%
% tfce_estimate_stat warns when the permutation null looks too narrow, by counting
% how many uncorrected p-values land in the upper tail. A check like that is only
% as good as its threshold and its direction, and this one had both wrong.
%
% WHICH UNIFORM. nPt is one-sided and conditioned on the sign of the observed
% effect: on mask_P it counts only the permutations with t >= t0, and t0 is
% positive there. Under the null the permuted statistic is symmetric about zero, so
% that probability can never exceed about a half. nPt is therefore uniform on
% (0, 0.5], NOT on (0, 1] -- an element carrying a large positive statistic cannot
% also be unusually SMALL. An earlier version compared nPt against 0.95 and expected
% ~5% to exceed it. None ever can, and none ever did.
%
% WHICH DIRECTION. A null that is too narrow makes every observed statistic look
% more extreme than it is, so the p-values are pushed DOWN and the upper tail
% EMPTIES. p_hi therefore FALLS below 5% when the null is too narrow, and rises
% above it when the null is too wide. The earlier version had the two branches the
% other way round.
%
% The checks below establish all of it: that the doubled p-values are uniform under
% a valid permutation, that the undoubled ones can never reach the old threshold,
% and that p_hi moves in the direction the warning now assumes.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Calibration of the permutation null');

rng(4);
n_vox  = 200000;
n_perm = 2000;

% ---------------------------------------------------------------------
% p_hi, computed exactly as tfce_estimate_stat computes it
%
% ratio is how much wider the TRUE distribution of the statistic is than the one
% the permutation scheme produces. 1 is a valid permutation; above 1 the
% permutation null is too narrow, which is the failure that matters.
% ---------------------------------------------------------------------
  function p_hi = calibration(ratio)
    t0 = ratio * randn(1, n_vox);            % the observed statistic
    tp = randn(n_perm, n_vox);               % what the permutations produce

    % one-sided, conditioned on the sign of the observed effect -- the same
    % construction as tperm/mask_P/mask_N in tfce_estimate_stat
    sgn = sign(t0); sgn(sgn == 0) = 1;
    cnt = sum(bsxfun(@times, tp, sgn) >= abs(t0), 1);
    nPt = cnt/n_perm;

    p_hi = mean(2*abs(nPt) > 0.95);          % the doubling is the fix
  end

% ---------------------------------------------------------------------
% under a valid permutation the doubled p-values are uniform
% ---------------------------------------------------------------------
p_valid = calibration(1.0);

val_util('result','a valid permutation gives ~5% in the upper tail', ...
  abs(p_valid - 0.05) < 0.01, ...
  sprintf('p_hi = %.3f (5%% expected)', p_valid));

% ---------------------------------------------------------------------
% the undoubled p-value can NEVER reach the threshold the old check used
% ---------------------------------------------------------------------
t0  = randn(1, n_vox);
tp  = randn(n_perm, n_vox);
sgn = sign(t0); sgn(sgn == 0) = 1;
nPt = sum(bsxfun(@times, tp, sgn) >= abs(t0), 1)/n_perm;

val_util('result','the p-value is one-sided, so it is bounded near 0.5', ...
  max(abs(nPt)) < 0.6, ...
  sprintf('largest uncorrected p = %.2f, not ~1', max(abs(nPt))));

val_util('result','the old threshold of 0.95 was unreachable', ...
  mean(abs(nPt) > 0.95) == 0, ...
  sprintf('%.1f%% of p exceed 0.95, so the check never fired', ...
          100*mean(abs(nPt) > 0.95)));

% ---------------------------------------------------------------------
% and p_hi moves the way the warning now assumes
% ---------------------------------------------------------------------
p_narrow = calibration(2.0);      % permutation null half as wide as it should be
p_wide   = calibration(0.5);      % twice as wide

fprintf('\n  %-34s %s\n','permutation null','p_hi');
fprintf('  %s\n', repmat('-',1,46));
fprintf('  %-34s %.3f\n','2x too narrow (anti-conservative)', p_narrow);
fprintf('  %-34s %.3f\n','valid',                             p_valid);
fprintf('  %-34s %.3f\n','2x too wide (conservative)',        p_wide);

val_util('result','a too-narrow null EMPTIES the upper tail', ...
  p_narrow < p_valid, ...
  sprintf('p_hi %.3f < %.3f, so the warning must fire on a LOW p_hi', ...
          p_narrow, p_valid));

val_util('result','a too-wide null FILLS it', ...
  p_wide > p_valid, ...
  sprintf('p_hi %.3f > %.3f', p_wide, p_valid));

% and the thresholds the code uses must actually separate the three cases
val_util('result','the 3% threshold fires on a too-narrow null', ...
  p_narrow < 0.03, sprintf('p_hi = %.3f', p_narrow));

val_util('result','and leaves a valid one alone', ...
  p_valid >= 0.03 && p_valid <= 0.08, sprintf('p_hi = %.3f', p_valid));

val_util('summary');

end
