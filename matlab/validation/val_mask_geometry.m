function val_mask_geometry
% Is mask/data geometry matching robust to harmless affine roundoff?
%
% ______________________________________________________________________
%
% Christian Gaser
% ______________________________________________________________________

val_util('header','Mask geometry check');
val_util('extract','geometry_mismatch');

V1.dim = [169 205 169];
V1.mat = [-1.5 0 0 126.75; 0 1.5 0 -145.5; 0 0 1.5 -72.75; 0 0 0 1];

V2 = V1;
val_util('result','identical geometry passes', ...
  ~geometry_mismatch(V1, V2, 1:3));

V2.mat(1,4) = V2.mat(1,4) + 5e-5;
val_util('result','small affine roundoff passes', ...
  ~geometry_mismatch(V1, V2, 1:3), ...
  'max |delta| = 5e-5 mm');

V2 = V1;
V2.mat(1,1) = V2.mat(1,1) + 5e-7;
val_util('result','small scale roundoff passes', ...
  ~geometry_mismatch(V1, V2, 1:3), ...
  'max corner displacement < 1e-4 mm');

V2 = V1;
V2.mat(1,4) = V2.mat(1,4) + 2e-4;
val_util('result','larger affine mismatch fails', ...
  geometry_mismatch(V1, V2, 1:3), ...
  'max |delta| = 2e-4 mm');

V2 = V1;
V2.mat(1,1) = V2.mat(1,1) + 1e-4;
val_util('result','larger scale mismatch fails', ...
  geometry_mismatch(V1, V2, 1:3), ...
  'max corner displacement > 1e-4 mm');

V2 = V1;
V2.dim(3) = V2.dim(3) + 1;
val_util('result','dimension mismatch fails', ...
  geometry_mismatch(V1, V2, 1:3));

M1.dim = 327684;
M1.mat = eye(4);
M2 = M1;
M2.mat(1,4) = M2.mat(1,4) + 5e-5;
val_util('result','mesh roundoff passes', ...
  ~geometry_mismatch(M1, M2, 1), ...
  'max corner displacement < 1e-4 mm');

M2 = M1;
M2.mat(1,4) = M2.mat(1,4) + 2e-4;
val_util('result','mesh mismatch fails', ...
  geometry_mismatch(M1, M2, 1), ...
  'max corner displacement > 1e-4 mm');

val_util('summary');
end
