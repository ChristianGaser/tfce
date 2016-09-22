tfce = tfceMex(t, deltaT, tbss, show_number_of_processors)
% FORMAT tfce = tfceMex(t, deltaT, show_number_of_processors)

% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling tfceMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O tfceMex.c
cd(p_path);

tfce = tfceMex(t, deltaT, tbss, show_number_of_processors);

return
