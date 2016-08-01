function make

disp('Compilerflag aendern in /TC')
% check for windows systems
if strcmp(mexext,'mexw32') || strcmp(mexext,'mexw64')
    obj = '.obj';
    WIN32 = ' -DSPM_WIN32 ';
else
    obj = '.o';
    WIN32 = '';
end

if strcmp(mexext,'mexmaci64')
   mexflag=' -O -largeArrayDims -Dchar16_t=UINT16_T CFLAGS=''$CFLAGS -pthread -Wall -ansi -pedantic -Wextra'' ';
else
   mexflag=' -O -largeArrayDims CFLAGS=''$CFLAGS -pthread -Wall -ansi -pedantic -Wextra'' ';
end

eval(['mex ' mexflag ' -O tfceMex_pthread.c'])

try
    tfceMex_pthread(rand(10,10,10),0.01,0.5,2);
    disp('Compilation of tfceMex successful')
catch
    disp('Compilation of tfceMex not successful')
end