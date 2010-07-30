function cg_tfce_estimate(SPM, Ic, xCon, n_perm0, vFWHM)

Pmat = spm_select(1,'SPM.mat','Select SPM.mat');
load(Pmat)

[Ic,xCon] = spm_conman(SPM,'T',Inf,...
        '  Select contrasts...',' for TFCE',1);

n_perm = 2500;
n_perm = spm_input('How many permutations? ',1,'r',n_perm,1,[10 n_perm]);
n_perm_break = spm_input('Stop point if no suprathreshold clusters are found? ',1,'r',100,1,[100 n_perm]);
vFWHM  = spm_input('Variance smoothing in FWHM (for low DFs) ','+1','e',0);

for k=1:length(Ic)
  cg_tfce_estimate(SPM, Ic(k), xCon, n_perm, vFWHM, n_perm_break)
end