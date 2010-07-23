function make

mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_CHAR 
movefile('spm_vol_utils.o', ['utils_uchar.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_SHORT 
movefile('spm_vol_utils.o', ['utils_short.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_INT 
movefile('spm_vol_utils.o', ['utils_int.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_CHAR 
movefile('spm_vol_utils.o', ['utils_schar.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT 
movefile('spm_vol_utils.o', ['utils_ushort.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_INT 
movefile('spm_vol_utils.o', ['utils_uint.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_FLOAT 
movefile('spm_vol_utils.o', ['utils_float.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_DOUBLE 
movefile('spm_vol_utils.o', ['utils_double.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP 
movefile('spm_vol_utils.o', ['utils_short_s.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP 
movefile('spm_vol_utils.o', ['utils_int_s.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP 
movefile('spm_vol_utils.o', ['utils_ushort_s.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP 
movefile('spm_vol_utils.o', ['utils_uint_s.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP 
movefile('spm_vol_utils.o', ['utils_float_s.' mexext '.o'],'f');
mex -O -c spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP 
movefile('spm_vol_utils.o', ['utils_double_s.' mexext '.o'],'f');
mex -O -c spm_make_lookup.c 
movefile('spm_make_lookup.o', ['spm_make_lookup.' mexext '.o'],'f');
mex -O -c spm_getdata.c 
movefile('spm_getdata.o', ['spm_getdata.' mexext '.o'],'f');
mex -O -c spm_vol_access.c 
movefile('spm_vol_access.o', ['spm_vol_access.' mexext '.o'],'f');
mex -O -c spm_mapping.c 
movefile('spm_mapping.o', ['spm_mapping.' mexext '.o'],'f');
eval(['mex -O cg_glm_get_Beta_ResSS.c utils_uchar.' mexext '.o utils_short.' mexext '.o utils_int.' mexext '.o utils_schar.' mexext '.o utils_ushort.' mexext '.o utils_uint.' mexext '.o utils_float.' mexext '.o utils_double.' mexext '.o utils_short_s.' mexext '.o utils_int_s.' mexext '.o utils_ushort_s.' mexext '.o utils_uint_s.' mexext '.o utils_float_s.' mexext '.o utils_double_s.' mexext '.o spm_make_lookup.' mexext '.o spm_getdata.' mexext '.o spm_vol_access.' mexext '.o spm_mapping.' mexext '.o']);

mex -O tfceMex.c