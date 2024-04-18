P = mfilename('fullpath');
[p,f,e]=spm_fileparts(P);
oldpwd = p;
cd(p);
mex spm_hist2_z_exp_polyval2.c
mex spm_hist2_masking.c
mex acid_c_dti_to_fa_HBM2010.c
%mex DTtoFA_2d.c
mex acid_c_dti_to_ev_ew.c


% POAS
cd ..
pPOAS = ['Preprocessing' filesep 'msPOAS'];
cd(pPOAS)
mex acid_c_poas_adsmse3ms.c
mex acid_c_poas_ghfse3i.c
mex acid_c_poas_ipolsp.c
mex acid_c_poas_linterpol.c
mex acid_c_poas_lkfse3i.c
mex acid_c_poas_lkfulls0.c
cd ..
cd ..
pHySCO = ['Preprocessing' filesep 'HySCO'];
cd(pHySCO)
mex acid_hysco_getPartialBMexC.cpp
mex acid_hysco_n2ccScalarMexC.cpp
cd('FAIRkernel')
mex linearInterMexC.cpp   
mex splineInterMexC.cpp
cd ..
cd ..
cd ..
% acid_hysco_EPImake('all')