tfce = tfceMex_pthread(t, dh, E, H, calc_neg, single_threaded)
% FORMAT tfce = tfceMex(t, dh, E, H, calc_neg, single_threaded)
% t         - T map 
% dh        - steps size (e.g. dh = max(abs(t))/100)
% E         - TFCE parameter
% H         - TFCE parameter
% calc_neg  - also calc neg. TFCE values (default)
% single_threaded - use single thread only
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling tfceMex_pthread.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O tfceMex_pthread.c
cd(p_path);

tfce = tfceMex_pthread(t, dh, E, H, calc_neg, single_threaded);

return
