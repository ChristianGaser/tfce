function compile

if ispc
  disp('Compilerflag aendern in /TC')
end

if strcmp(mexext,'mexmaci64') && verLessThan('matlab','9.2')
  warning('WARNING: Matlab version should be at least R2017a for compilation under Mac.');
end

mexflag=' -O -largeArrayDims COPTIMFLAGS=''-O3 -fwrapv -DNDEBUG'' CFLAGS=''$CFLAGS -pthread -Wall -ansi -pedantic -Wextra'' ';

eval(['mex ' mexflag ' tfceMex_pthread.c'])

try
    tfceMex_pthread(rand(10,10,10),0.01,0.5,2);
    disp('Compilation of tfceMex successful')
catch
    disp('Compilation of tfceMex not successful')
end