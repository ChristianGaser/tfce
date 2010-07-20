tfce = tfceMex(t, n_steps)
% FORMAT tfce = tfceMex(t, n_steps)

% Christian Gaser
% $Id: tfceMex.m 224 2009-12-02 23:39:15Z gaser $

rev = '$Rev: 224 $';

disp('Compiling tfceMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O tfceMex.c
cd(p_path);

tfce = tfceMex(t, n_steps)

return
