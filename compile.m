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
   mexflag=' -O -largeArrayDims -Dchar16_t=UINT16_T CFLAGS=''$CFLAGS -Wall -ansi -pedantic -Wextra'' CPPLAGS=''$CPPFLAGS -Wall -ansi -pedantic -Wextra'' ';
else
    mexflag=' -O -largeArrayDims ';
end

if 0
mex  -O -c spm_vol_utils.c -DSPM_UNSIGNED_CHAR 
movefile(['spm_vol_utils' obj], ['utils_uchar.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_SIGNED_SHORT 
movefile(['spm_vol_utils' obj], ['utils_short.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_SIGNED_INT 
movefile(['spm_vol_utils' obj], ['utils_int.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_SIGNED_CHAR 
movefile(['spm_vol_utils' obj], ['utils_schar.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT 
movefile(['spm_vol_utils' obj], ['utils_ushort.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_UNSIGNED_INT 
movefile(['spm_vol_utils' obj], ['utils_uint.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_FLOAT 
movefile(['spm_vol_utils' obj], ['utils_float.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_DOUBLE 
movefile(['spm_vol_utils' obj], ['utils_double.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_short_s.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_int_s.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_ushort_s.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_uint_s.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_float_s.' mexext obj],'f');
mex  -I/usr/include -O -c spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP 
movefile(['spm_vol_utils' obj], ['utils_double_s.' mexext obj],'f');
mex  -I/usr/include -O -c spm_make_lookup.c 
movefile(['spm_make_lookup' obj], ['spm_make_lookup.' mexext obj],'f');
mex  -I/usr/include -O -c spm_getdata.c 
movefile(['spm_getdata' obj], ['spm_getdata.' mexext obj],'f');
eval(['mex -I/usr/include ' WIN32 ' -O -c spm_vol_access.c']) 
movefile(['spm_vol_access' obj], ['spm_vol_access.' mexext obj],'f');

eval(['mex ' WIN32 mexflag ' -O -c spm_mapping.c'])
movefile(['spm_mapping' obj], ['spm_mapping.' mexext obj],'f');
end

try % try OpenMP support
    if strcmp(mexext,'mexmaci64')
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-U_OPENMP -m64 -fPIC -O3'' -I/usr/include ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
        movefile(['cg_glm_get_Beta_ResSS.' mexext], ['cg_glm_get_Beta_ResSS_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-U_OPENMP -m64 -fPIC -O3' -I/usr/include -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-fopenmp -m64 -fPIC -O3' -O /usr/local/lib/x86_64/libgomp.a tfceMex.c         
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-fopenmp -m64 -fPIC -O3'' -I/usr/include ' WIN32 ' -O /usr/local/lib/x86_64/libgomp.a cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
    elseif strcmp(mexext,'mexmaci')
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-U_OPENMP -m32 -fPIC -O3'' -I/usr/include ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
        movefile(['cg_glm_get_Beta_ResSS.' mexext], ['cg_glm_get_Beta_ResSS_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-U_OPENMP -m32 -fPIC -O3' -I/usr/include -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/local/lib/x86/libgomp.a tfceMex.c 
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-fopenmp -m32 -fPIC -O3'' -I/usr/include ' WIN32 ' -O /usr/local/lib/x86/libgomp.a cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
    elseif strcmp(mexext,'mexa64')
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-U_OPENMP -m64 -fPIC -O3'' -I/usr/include ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
        movefile(['cg_glm_get_Beta_ResSS.' mexext], ['cg_glm_get_Beta_ResSS_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-U_OPENMP -m64 -fPIC -O3' -I/usr/include -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-fopenmp -m64 -fPIC -O3' -O -lgomp tfceMex.c
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-fopenmp -m64 -fPIC -O3'' -I/usr/include ' WIN32 ' -lgomp -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
    elseif strcmp(mexext,'mexglx')
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-U_OPENMP -m32 -fPIC -O3'' -I/usr/include ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
        movefile(['cg_glm_get_Beta_ResSS.' mexext], ['cg_glm_get_Beta_ResSS_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-U_OPENMP -m32 -fPIC -O3' -I/usr/include -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.4' CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/lib/i386-linux-gnu/gcc/i686-linux-gnu/4.4/libgomp.a tfceMex.c 
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-fopenmp -m32 -fPIC -O3'' -I/usr/include ' WIN32 ' -O /usr/lib/i386-linux-gnu/gcc/i686-linux-gnu/4.4/libgomp.a cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
    elseif strcmp(mexext,'mexw64')
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-U_OPENMP -m64 -fPIC -O3'' ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
        movefile(['cg_glm_get_Beta_ResSS.' mexext], ['cg_glm_get_Beta_ResSS_noopenmp.' mexext],'f');
        mex CFLAGS='-U_OPENMP -m64 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp -m64 -fPIC -O3' -O tfceMex.c
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-fopenmp -m64 -fPIC -O3'' ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
    elseif strcmp(mexext,'mexw32')
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-U_OPENMP -m32 -fPIC -O3'' ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
        movefile(['cg_glm_get_Beta_ResSS.' mexext], ['cg_glm_get_Beta_ResSS_noopenmp.' mexext],'f');
        mex CFLAGS='-U_OPENMP -m32 -fPIC -O3' -O tfceMex.c
        movefile(['tfceMex.' mexext], ['tfceMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp -m32 -fPIC -O3' -O tfceMex.c
        eval(['mex CC=''gcc-4.4'' CFLAGS=''-fopenmp -m32 -fPIC -O3'' ' WIN32 ' -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext obj ' utils_short.' mexext obj ...
            ' utils_int.' mexext obj ' utils_schar.' mexext obj ' utils_ushort.' mexext obj ' utils_uint.' mexext obj ...
            ' utils_float.' mexext obj ' utils_double.' mexext obj ' utils_short_s.' mexext obj ' utils_int_s.' mexext obj ...
            ' utils_ushort_s.' mexext obj ' utils_uint_s.' mexext obj ' utils_float_s.' mexext obj ' utils_double_s.' mexext obj ...
            ' spm_make_lookup.' mexext obj ' spm_getdata.' mexext obj ' spm_vol_access.' mexext obj ' spm_mapping.' mexext obj]);
    end
    disp('Compiling tfceMex with OpenMP')
catch 
    disp('Compiling tfceMex without OpenMP')
    mex CC='gcc-4.4' CFLAGS='-fPIC -O3' -O tfceMex.c 
end

try
    tfceMex(rand(10,10,10),5);
    disp('Compilation of tfceMex successful')
catch
    disp('Compilation of tfceMex not successful')
end