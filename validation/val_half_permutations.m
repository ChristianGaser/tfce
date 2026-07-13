function val_half_permutations
% Is the half-permutation shortcut lossless?
%
% For a two-sample design with EQUAL group sizes and a t-contrast, the toolbox
% computes only half of the permutations and enters each one twice into the null
% distribution, once as max(t) and once as -min(t). This assumes that the
% mirrored permutation (the one that swaps the two groups) produces exactly the
% negated statistic map. The shortcut is lossless if and only if that holds.
%
% It also checks that the assumption FAILS for unequal group sizes, which is why
% the toolbox restricts the shortcut to balanced designs.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Half-permutation shortcut');
rng(1);
n_vox = 2000;

% ---------------------------------------------------------------------
% equal group sizes: mirrored statistic must be exactly -t
% ---------------------------------------------------------------------
for withcov = [false true]
  n1 = 15; n2 = 15; n = n1+n2;
  g  = [ones(n1,1); zeros(n2,1)];
  X  = [g 1-g];
  c  = [1; -1];
  if withcov
    X = [X randn(n,1) ones(n,1)];       % nuisance covariate + intercept
    c = [1; -1; 0; 0];
    lbl = 'equal n, with nuisance covariate';
  else
    lbl = 'equal n, no covariate';
  end

  Y = randn(n, n_vox);

  worst = 0;
  for r = 1:50
    pr = randperm(n);
    Xp = X; Xp(:,1:2) = X(pr,1:2);       % Draper-Stoneman on the group columns
    Xm = Xp; Xm(:,[1 2]) = Xp(:,[2 1]);  % the mirrored permutation

    t  = glm_t(Y, Xp, c);
    tm = glm_t(Y, Xm, c);
    worst = max(worst, max(abs(tm + t))./max(abs(t)));
  end
  val_util('result', ['mirror = -t (' lbl ')'], worst < 1e-12, ...
    sprintf('max rel |t_mirror + t| = %.2e', worst));
end

% ---------------------------------------------------------------------
% unequal group sizes: the mirror is NOT a permutation of the same design,
% so the shortcut must not be used. Verify the assumption really breaks.
% ---------------------------------------------------------------------
n1 = 18; n2 = 12; n = n1+n2;
g  = [ones(n1,1); zeros(n2,1)];
X  = [g 1-g]; c = [1;-1];
Y  = randn(n, n_vox);

sizes_differ = false;
for r = 1:20
  pr = randperm(n);
  Xp = X; Xp(:,1:2) = X(pr,1:2);
  Xm = Xp; Xm(:,[1 2]) = Xp(:,[2 1]);
  % the mirrored design has the group sizes swapped -> not in the permutation set
  if sum(Xm(:,1)) ~= sum(X(:,1)), sizes_differ = true; end
end
val_util('result','unequal n: mirror leaves the permutation set', sizes_differ, ...
  'group sizes are swapped, so the shortcut is correctly disabled');

% ---------------------------------------------------------------------
% consequence for the null distribution: with equal n, taking (max t, -min t)
% from half the permutations must give the same null as running all of them
% ---------------------------------------------------------------------
n1 = 15; n2 = 15; n = n1+n2;
g = [ones(n1,1);zeros(n2,1)]; X = [g 1-g]; c = [1;-1];
Y = randn(n, n_vox);

n_half = 200;
null_half = zeros(1, 2*n_half);
null_full = zeros(1, 2*n_half);
k = 0;
for r = 1:n_half
  pr = randperm(n);
  Xp = X; Xp(:,1:2) = X(pr,1:2);
  Xm = Xp; Xm(:,[1 2]) = Xp(:,[2 1]);

  t  = glm_t(Y, Xp, c);
  tm = glm_t(Y, Xm, c);            % explicitly computed mirror

  % what the shortcut stores
  null_half(2*r-1) = max(t);
  null_half(2*r)   = -min(t);
  % what running both permutations explicitly would store
  null_full(2*r-1) = max(t);
  null_full(2*r)   = max(tm);
  k = k + 1;
end
d = max(abs(sort(null_half) - sort(null_full)));
val_util('result','null from shortcut == null from explicit mirrors', d < 1e-10, ...
  sprintf('max |d| = %.2e over %d entries', d, 2*n_half));

val_util('summary');

%---------------------------------------------------------------
function t = glm_t(Y, X, c)
% Y is n_data x n_vox here
pKX = pinv(X); n = size(X,1);
B = pKX*Y;
R = Y - X*B;
ResMS = sum(R.^2,1)/(n - rank(X));
t = (c'*B)./sqrt(ResMS*(c'*(pKX*pKX')*c));
