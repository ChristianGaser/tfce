function make

% check for windows systems
if strcmp(mexext,'mexw32') || strcmp(mexext,'mexw64')
    obj = '.obj';
    WIN32 = '-DSPM_WIN32 ';
else
    obj = '.o';
    WIN32 = '';
end

mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_CHAR 
movefile(['spm_vol_utils' obj], ['utils_uchar.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_SHORT 
movefile(['spm_vol_utils' obj], ['utils_short.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_INT 
movefile(['spm_vol_utils' obj], ['utils_int.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_CHAR 
movefile(['spm_vol_utils' obj], ['utils_schar.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT 
movefile(['spm_vol_utils' obj], ['utils_ushort.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_INT 
movefile(['spm_vol_utils' obj], ['utils_uint.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_FLOAT 
movefile(['spm_vol_utils' obj], ['utils_float.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_DOUBLE 
movefile(['spm_vol_utils' obj], ['utils_double.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_short_s.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_int_s.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_ushort_s.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_uint_s.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_float_s.' mexext obj],'f');
mex -O -c spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_double_s.' mexext obj],'f');
mex -O -c spm_make_lookup.c 
movefile(['spm_make_lookup' obj], ['spm_make_lookup.' mexext obj],'f');
mex -O -c spm_getdata.c 
movefile(['spm_getdata' obj], ['spm_getdata.' mexext obj],'f');
eval(['mex ' WIN32 ' -O -c spm_vol_access.c']) 
movefile(['spm_vol_access' obj], ['spm_vol_access.' mexext obj],'f');
eval(['mex ' WIN32 ' -O -c spm_mapping.c']) 
movefile(['spm_mapping' obj], ['spm_mapping.' mexext obj],'f');
eval(['mex ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
    ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
    ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
    ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
    ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);

try % try OpenMP support
    if strcmp(mexext,'mexmaci64')
        mex CC='gcc-4.2' CFLAGS='-m64 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.2' CFLAGS='-fopenmp -m64 -fPIC -O3' -O /usr/local/lib/x86_64/libgomp.a tfceMex.c 
    elseif strcmp(mexext,'mexmaci')
        mex CC='gcc-4.2' CFLAGS='-m32 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.2' CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/local/lib/x86/libgomp.a tfceMex.c 
    elseif strcmp(mexext,'mexa64')
        mex CFLAGS='-m64 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp -m64 -fPIC -O3' -O -lgomp tfceMex.c
    elseif strcmp(mexext,'mexglx')
        mex CFLAGS='-m32 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/lib/gcc/i486-linux-gnu/4.4/libgomp.a tfceMex.c 
    elseif strcmp(mexext,'mexw64')
        mex CFLAGS='-m64 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp m64 -fPIC -O3' -O tfceMex.c
    elseif strcmp(mexext,'mexw32')
        mex CFLAGS='-m32 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp m32 -fPIC -O3' -O tfceMex.c
    end
    disp('Compiling tfceMex with OpenMP')
catch 
    disp('Compiling tfceMex without OpenMP')
    mex CFLAGS='-fPIC -O3' -O tfceMex.c 
end

try
    tfceMex(rand(10,10,10),5);
    disp('Compilation of tfceMex successful')
catch
    disp('Compilation of tfceMex not successful')
end