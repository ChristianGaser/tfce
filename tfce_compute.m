function tfce = tfce_compute(t, E, H, calc_neg, faces)
% Apply the TFCE transform, either exactly or with the dh-stepping approximation
% FORMAT tfce = tfce_compute(t, E, H, calc_neg, opt)
% t         - T map (3D volume, or Nx1 vector for surface data)
% E         - TFCE parameter for extent
% H         - TFCE parameter for height
% calc_neg  - also calc neg. TFCE values
% faces     - surface faces; empty for volume data
%
% The dh-stepping variant approximates the TFCE integral on a grid of n_steps
% levels, with dh derived from the map maximum. The max-tree variant evaluates
% the integral exactly and has no step size, so n_steps and singlethreaded are
% ignored.
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if isempty(faces)
  tfce = tfceMex_maxtree(t, E, H, calc_neg);
else
  tfce = tfceMex_maxtree(t, E, H, calc_neg, faces);
end
