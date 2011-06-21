tfce = tfceMex(t, n_steps, show_number_of_processors)
% FORMAT tfce = tfceMex(t, n_steps, show_number_of_processors)

% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling tfceMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O tfceMex.c
cd(p_path);

tfce = tfceMex(t, n_steps, show_number_of_processors);

return
