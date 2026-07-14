function val_tfce_exactness
% Is the max-tree TFCE the exact value of the TFCE integral?
%
% The mex-file is compared against an INDEPENDENT reference implementation
% written in plain MATLAB, which approximates
%
%   TFCE(v) = int_0^{t_v} e_v(h)^E h^H dh
%
% by stepping through n thresholds and labelling connected components with
% bwlabeln / the mesh adjacency at each step. As the step size goes to zero the
% reference must converge onto the max-tree result, and it must do so at first
% order (the error halving whenever n is doubled). This establishes that the
% max-tree computes the limit of the conventional approximation rather than a
% different quantity.
%
% ______________________________________________________________________
% Christian Gaser
% ______________________________________________________________________

val_util('header','Exactness of the max-tree TFCE');

rng(4);
E = 0.5; H = 2.0;

% ---------------------------------------------------------------------
% volume
% ---------------------------------------------------------------------
dim = [40 40 40];
t = smooth3(randn(dim),'gaussian',7,1.8);
t = t/std(t(:));
[x,y,z] = ndgrid(linspace(-1,1,dim(1)),linspace(-1,1,dim(2)),linspace(-1,1,dim(3)));
t = t + 2.5*exp(-((x-0.2).^2+(y-0.1).^2+(z+0.1).^2)/0.02);
t = t - 2.0*exp(-((x+0.3).^2+(y+0.25).^2+(z-0.2).^2)/0.02);

exact = tfceMex_maxtree(t, E, H, 1);

fprintf('\n  volume (%dx%dx%d), reference = plain MATLAB dh-stepping\n', dim);
fprintf('  %-10s %-14s\n','n_steps','max rel err');
ns = [50 100 200 400]; errs = zeros(size(ns));
for k = 1:numel(ns)
  ref = tfce_reference_vol(t, E, H, ns(k));
  errs(k) = max(abs(ref(:)-exact(:)))/max(abs(exact(:)));
  fprintf('  %-10d %-14.3e\n', ns(k), errs(k));
end
% First-order convergence means err ~ C/n, i.e. slope of log(err) vs log(n) = -1.
% This tests the trend rather than each individual ratio, which is noisy because
% the error is a maximum over elements.
pf = polyfit(log(ns), log(errs), 1); slope = pf(1);
val_util('result','volume: dh-stepping converges onto the max-tree', ...
  slope > -1.2 && slope < -0.8, ...
  sprintf('log-log slope = %.2f (first order = -1)', slope));

% ---------------------------------------------------------------------
% surface
% ---------------------------------------------------------------------
if exist('spm_mesh_adjacency','file') && exist('gifti','file')
  g = gifti(fullfile(spm('dir'),'canonical','cortex_5124.surf.gii'));
  faces = double(g.faces); nv = size(g.vertices,1);
  A = spm_mesh_adjacency(faces) + speye(nv);
  D = spdiags(1./sum(A,2),0,nv,nv)*A;
  s = randn(nv,1); for i=1:10, s = D*s; end
  s = s/std(s); s(1:round(nv/20)) = s(1:round(nv/20)) + 3;

  exact_s = tfceMex_maxtree(s, E, H, 1, faces);

  fprintf('\n  surface (%d vertices)\n', nv);
  fprintf('  %-10s %-14s\n','n_steps','max rel err');
  ns = [50 100 200 400]; errs = zeros(size(ns));
  for k = 1:numel(ns)
    ref = tfce_reference_surf(s, spm_mesh_adjacency(faces), E, H, ns(k));
    errs(k) = max(abs(ref(:)-exact_s(:)))/max(abs(exact_s(:)));
    fprintf('  %-10d %-14.3e\n', ns(k), errs(k));
  end
  pf = polyfit(log(ns), log(errs), 1); slope = pf(1);
  val_util('result','surface: dh-stepping converges onto the max-tree', ...
    slope > -1.2 && slope < -0.8, ...
    sprintf('log-log slope = %.2f (first order = -1)', slope));

  % ------------------------------------------------------------------
  % the faces must be read in the class they actually arrive in
  %
  % GIFTI stores triangles as int32, so SPM.xVol.G.faces is an int32 array and
  % that is exactly what tfce_estimate_stat passes. Reading it as double -- which
  % is what mxGetPr does, whatever the array really holds -- reinterprets the
  % integer bit patterns as floating point and yields vertex indices that are
  % noise. Those were then written past the ends of the adjacency arrays: a heap
  % corruption on every single surface analysis.
  %
  % Note the cast on line 58 above. Testing only with double faces is precisely
  % how this was missed, so the class is now part of the test.
  % ------------------------------------------------------------------
  faces_i32 = int32(g.faces);
  val_util('result','surface: faces really are int32 as SPM passes them', ...
    isa(g.faces,'int32'), sprintf('class(g.faces) = %s', class(g.faces)));

  s_i32 = tfceMex_maxtree(s, E, H, 1, faces_i32);
  val_util('result','surface: int32 faces give the same map as double faces', ...
    isequal(s_i32, exact_s), 'bit-identical');

  b_i32 = tfceMex_maxtree_batch([s s], E, H, 1, faces_i32, 2);
  val_util('result','surface: the batched transform accepts int32 faces too', ...
    isequal(b_i32(:,1), exact_s(:)), 'and matches the single-map transform');

  % an index naming a vertex that does not exist must be refused, not written
  % past the end of the adjacency arrays
  bad = faces_i32; bad(1) = int32(nv + 5);
  refused = false;
  try
    tfceMex_maxtree(s, E, H, 1, bad);
  catch
    refused = true;
  end
  val_util('result','surface: an out-of-range vertex index is refused', ...
    refused, 'rather than corrupting the heap');
else
  fprintf('\n  (SPM not on the path, skipping the surface check)\n');
end

% ---------------------------------------------------------------------
% the batched transform must be identical to the single-map one
% ---------------------------------------------------------------------
B = 8; N = prod(dim);
T = zeros(N,B);
for b = 1:B
  m = smooth3(randn(dim),'gaussian',7,1.8);
  T(:,b) = m(:);
end
S = tfceMex_maxtree_batch(T, E, H, 1, dim);
d = 0;
for b = 1:B
  s1 = tfceMex_maxtree(reshape(T(:,b),dim), E, H, 1);
  d = max(d, max(abs(s1(:)-S(:,b))));
end
val_util('result','batched transform == single-map transform', d == 0, ...
  sprintf('max abs diff = %.1e over %d maps', d, B));

% ---------------------------------------------------------------------
% the permutation loop hands a whole block of permutations to the batched
% transform under ONE calc_neg flag, whereas it used to decide that per map. That
% is only the same thing if asking for negative TFCE values on a map that has no
% negative values leaves it unchanged -- otherwise a block that mixes signed and
% unsigned maps (an F-statistic, a t-statistic that happens to be all-positive)
% would come out differently than it did permutation by permutation.
% ---------------------------------------------------------------------
T = zeros(N,B);
for b = 1:B
  m = smooth3(randn(dim),'gaussian',7,1.8);
  if mod(b,2) == 0, m = abs(m); end     % half the block has no negative values
  T(:,b) = m(:);
end

% sequential, deciding calc_neg per map, exactly as the loop used to
seq = zeros(N,B);
for b = 1:B
  t = reshape(T(:,b),dim);
  s1 = tfceMex_maxtree(t, E, H, min(t(:)) < 0);
  seq(:,b) = s1(:);
end

% batched, with one flag for the whole block, as the loop does now
bat = tfceMex_maxtree_batch(T, E, H, any(T(:) < 0), dim);

val_util('result','one calc_neg flag for a mixed-sign block is the same thing', ...
  isequal(seq, bat), ...
  sprintf('%d of %d maps had no negative values', nnz(~any(T<0,1)), B));

val_util('summary');

%---------------------------------------------------------------
function tfce = tfce_reference_vol(t, E, H, n_steps)
% independent reference: step through thresholds, label components with bwlabeln
tfce = zeros(size(t));
dh = max(abs(t(:)))/n_steps;
for sgn = [1 -1]
  d = sgn*t;
  for i = 1:n_steps
    h = i*dh;
    [L, num] = bwlabeln(d >= h, 26);
    if num == 0, continue; end
    sz = accumarray(L(L>0), 1);
    idx = L > 0;
    tfce(idx) = tfce(idx) + sgn * (sz(L(idx)).^E) * (h^H) * dh;
  end
end

%---------------------------------------------------------------
function tfce = tfce_reference_surf(t, A, E, H, n_steps)
% independent reference for surfaces, using the mesh adjacency
t = t(:);
tfce = zeros(size(t));
dh = max(abs(t))/n_steps;
for sgn = [1 -1]
  d = sgn*t;
  for i = 1:n_steps
    h = i*dh;
    sup = find(d >= h);
    if isempty(sup), continue; end
    As = A(sup,sup);
    [nc, C] = graphconncomp_local(As);
    if nc == 0, continue; end
    sz = accumarray(C(:), 1);
    tfce(sup) = tfce(sup) + sgn * (sz(C(:)).^E) * (h^H) * dh;
  end
end

%---------------------------------------------------------------
function [nc, C] = graphconncomp_local(A)
% connected components of a sparse symmetric adjacency, via Dulmage-Mendelsohn
A = A + speye(size(A));
[p,~,r] = dmperm(A);
nc = numel(r) - 1;
C = zeros(size(A,1),1);
for i = 1:nc
  C(p(r(i):r(i+1)-1)) = i;
end
