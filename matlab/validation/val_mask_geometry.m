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
V2.mat(1,4) = V2.mat(1,4) + 2e-4;
val_util('result','larger affine mismatch fails', ...
  geometry_mismatch(V1, V2, 1:3), ...
  'max |delta| = 2e-4 mm');

V2 = V1;
V2.dim(3) = V2.dim(3) + 1;
val_util('result','dimension mismatch fails', ...
  geometry_mismatch(V1, V2, 1:3));

val_util('summary');
end
