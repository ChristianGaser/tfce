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

eval(['mex ' mexflag ' tfceMex_maxtree.c'])

try
    tfceMex_maxtree(rand(10,10,10),0.5,2,1);
    disp('Compilation of tfceMex_maxtree successful')
catch
    disp('Compilation of tfceMex_maxtree not successful')
end

eval(['mex ' mexflag ' tfceMex_maxtree_batch.c'])

try
    tfceMex_maxtree_batch(rand(1000,4),0.5,2,1,[10 10 10]);
    disp('Compilation of tfceMex_maxtree_batch successful')
catch
    disp('Compilation of tfceMex_maxtree_batch not successful')
end

eval(['mex ' mexflag ' tfceMex_resss.c'])

try
    tfceMex_resss(single(rand(100,10)),single(rand(100,2)),rand(10,2));
    disp('Compilation of tfceMex_resss successful')
catch
    disp('Compilation of tfceMex_resss not successful')
end