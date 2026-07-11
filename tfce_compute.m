function tfce = tfce_compute(t, E, H, calc_neg, opt)
% Apply the TFCE transform, either exactly or with the dh-stepping approximation
% FORMAT tfce = tfce_compute(t, E, H, calc_neg, opt)
% t         - T map (3D volume, or Nx1 vector for surface data)
% E         - TFCE parameter for extent
% H         - TFCE parameter for height
% calc_neg  - also calc neg. TFCE values
% opt       - options with fields
%             use_maxtree    - exact max-tree TFCE instead of dh-stepping
%             n_steps        - number of dh steps (dh-stepping only)
%             faces          - surface faces; empty for volume data
%             singlethreaded - dh-stepping volumes only
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

if opt.use_maxtree
  if isempty(opt.faces)
    tfce = tfceMex_maxtree(t, E, H, calc_neg);
  else
    tfce = tfceMex_maxtree(t, E, H, calc_neg, opt.faces);
  end
else
  dh = max(abs(t(:)))/opt.n_steps;
  if isempty(opt.faces)
    tfce = tfceMex_pthread(t, dh, E, H, calc_neg, opt.singlethreaded)*dh;
  else
    tfce = tfce_mesh(opt.faces, t, dh, E, H)*dh;
  end
end
