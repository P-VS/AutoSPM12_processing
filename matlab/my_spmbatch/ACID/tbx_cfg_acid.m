function acid = tbx_cfg_acid

fprintf('\n');
disp( 'Initialising ACID                                                                              ');
disp( '    _     ____   _   ____                                                                      ');
disp( '   /_\   / ___| | | |  _  \   A Comprehensive Toolbox for Image Processing and Modeling        ');
disp( '  / _ \ | |___  | | | |_| |   of Brain, Spinal Cord, and Ex Vivo Diffusion MRI Data            ');
disp( ' /_/ \_\ \____| |_| |____ /   ACID Toolbox -- http://www.diffusiontools.com                    ');
disp( '                                                                                               ');





% disp( 'Initialising ACID                                                                               ');
% disp( '    _      ____   _   ____                                                                      ');
% disp( '   / \    / ___| | | |  _  \                                                                    ');
% disp( '  / _ \  | |     | | | | | |   A Comprehensive Toolbox for Image Processing and Modeling        ');
% disp( ' / ___ \ | |___  | | | |_| |   of Brain, Spinal Cord, and Post-mortem Diffusion MRI Data        ');
% disp( '/_/   \_\ \____| |_| |____ /   ACID-Toolbox -- https://www.diffusiontools.com                   ');
% disp( '                                                                                                ');

%% Add paths
d = fileparts(mfilename('fullpath'));
if ~isdeployed
    
    addpath(fullfile(d));

    %% read image
    addpath(fullfile(d,'read_write_data'));

    %% TPM
    addpath(fullfile(d,'ACID_TPM'));

    %% pre-processing
    addpath(fullfile(d,'Preprocessing'));
    addpath(fullfile(d,'Utilities'));
    pu = [d filesep 'Utilities'];
    addpath(fullfile(pu,'image_manipulation'));
    addpath(fullfile(pu,'mask_generation'));
    addpath(fullfile(pu,'quality_assessment'));
    pd = [d filesep 'Preprocessing'];
    addpath(fullfile(pd,'ECMOCO'));
    addpath(fullfile(pd,'RBC'));
    addpath(fullfile(pd,'koays-inversion'));
    pk = genpath(fullfile(pd,'koays-inversion'));
    addpath(pk);
    addpath(fullfile(pd,'md-dmri'));
    pk = genpath(fullfile(pd,'md-dmri'));
    addpath(pk);
    pr = [pd filesep 'RBC'];
    addpath(fullfile(pr,'rician-bias-simulation'));
    %% hysco
    if exist(fullfile(pd,'HySCO'),'dir')
        addpath(fullfile(pd,'HySCO'));
        path_hysco = [pd filesep 'HySCO'];
        addpath(fullfile(path_hysco,'FAIRkernel'));
        dummy_hysco = 1;
    else
        dummy_hysco = 0;
    end

    %% msPOAS
    if exist(fullfile(pd,'msPOAS'),'dir')
        addpath(fullfile(pd,'msPOAS'));
    end

    %% DTI
    dummy_lpf = 0;

    addpath(fullfile(d,'cfiles'));
    addpath(fullfile(d,'TensorFit'));
    addpath(fullfile(d,'config'));
    ld = [d filesep 'config'];
    addpath(fullfile(ld,'local'));
    addpath(fullfile(d,'TensorFit','Gauss_Newton_Fit'));

    %% Kurtosis estimation
    Tensorfitpath = [d filesep 'TensorFit'];
    addpath(fullfile(Tensorfitpath,'Kurtosis'));
    addpath(fullfile(Tensorfitpath,['Kurtosis' filesep 'Elliptic_Integrals']));
    addpath(fullfile(Tensorfitpath,['Kurtosis' filesep 'quadraticprogramming']));
    addpath(fullfile(Tensorfitpath,'writeData'));

    %% biophys
    if exist(fullfile(d,'biophys'),'dir')
        addpath(fullfile(d,'biophys'));
    end

    %% NODDI-DTI
    if exist(fullfile(d,'biophys','noddi_dti'),'dir')
        addpath(fullfile(d,'biophys','noddi_dti'));
        dummy_noddidti =1;
    else
        dummy_noddidti =0;
    end

    %% WMTI-WATSON
    if exist(fullfile(d,'TensorFit','wmti'),'dir')
       addpath(fullfile(d,'TensorFit','wmti')); 
    end
    %% WMTI
    if exist(fullfile(d,'biophys','wmti_watson'),'dir')
       addpath(fullfile(d,'biophys','wmti_watson')); 
    end
    %% Tensor Fiber Density - REiser et al., NI, 2013
    if exist(fullfile(d,'Freiburg_DTItools'),'dir')
        df = [d filesep 'Freiburg_DTItools'];
        addpath(fullfile(d,'Freiburg_DTItools'));
        if exist(fullfile(df,'DTI_tools'),'dir')
            addpath(fullfile(df,'DTI_tools'));
        end
        if exist(fullfile(df,'impexp'),'dir')
            addpath(fullfile(df,'impexp'));
        end
    end

    % check whether multi-hysco-field option exist
    path_hysco = [pd filesep 'HySCO'];
    if exist(fullfile(path_hysco,'acid_hysco_write_multiB.m'),'file')
        dummy_hysco_apply_multi = true;
    else
        dummy_hysco_apply_multi = false;
    end
else
    %% Example_Batches
    addpath(fullfile(d,'Example_Batches'));

    %% configuration
    addpath(fullfile(d,'config'));
    ld = [d filesep 'config'];
    addpath(fullfile(ld,'local'));

    %% pre-processing
    pd = [d filesep 'Preprocessing'];
    if(exist(fullfile(pd,'HySCO'),'dir'))
        dummy_hysco = 1;
    else
        dummy_hysco = 0;
    end
    
    % check whether multi-hysco-field option exist
    path_hysco = [d filesep 'HySCO'];
    if exist(fullfile(path_hysco, 'HySCO_write_multiB.m'),'file')
        dummy_hysco_apply_multi = true;
    else
        dummy_hysco_apply_multi = false;
    end

    %% DTI
    dummy_lpf = 0;

    %% NODDI-DTI
    if exist(fullfile(pd,'biophys','noddidti'),'dir')
        dummy_noddidti = 1;
    else
        dummy_noddidti = 0;
    end
end

% configuration
addpath(fullfile(d,'config'));
ld = [d filesep 'config'];
addpath(fullfile(ld,'local'));
% config = tbx_main_acid_config;

% ACID_Test
addpath(fullfile(d,'ACID_Test'));
td = [d filesep 'ACID_Test'];
addpath(fullfile(td,'Batch_Files'));
addpath(fullfile(td,'Fit_Results_Of_Test'));
addpath(fullfile(td,'Ground_Truth_Data'));

%% LPF correction
if dummy_lpf
    % order of spherical harmonics
    Nl_lpf        = cfg_menu;
    Nl_lpf.tag    = 'Nl_lpf';
    Nl_lpf.name   = 'Order of spherical harmonics';
    Nl_lpf.help   = {'Choose the order of spherical harmonics. Note that currently only spherical harmonics up to the order of three are available.'};
    Nl_lpf.labels = {
        '1st order spherical harmonics'
        '2nd order spherical harmonics'
        '3rd order spherical harmonics'
        }';
    Nl_lpf.values = {0 1 2};
    Nl_lpf.val    = {2};

    % regularisation factor epsilon
    epsilon_lpf         = cfg_entry;
    epsilon_lpf.tag     = 'epsilon_lpf';
    epsilon_lpf.name    = 'regularisation factor';
    epsilon_lpf.help    = {'Regularisation factor. See Mohammadi et al., Neuroimage 2012 for details. Do not touch default value if you are not sure.'};
    epsilon_lpf.strtype = 'e';
    epsilon_lpf.num     = [1 1];
    epsilon_lpf.val     = {1/10};

    % percentage coverage of brain mask
    perc_lpf         = cfg_entry;
    perc_lpf.tag     = 'perc_lpf';
    perc_lpf.name    = 'brain mask parameter';
    perc_lpf.help    = {'Factor that depends on ratio between brain coverage and field-of-view. Less brain coverage of field-of-view means lower perc-value.'};
    perc_lpf.strtype = 'e';
    perc_lpf.num     = [1 1];
    perc_lpf.val     = {0.8};

    % b values
    b_vals_lpf         = cfg_entry;
    b_vals_lpf.tag     = 'b_vals_lpf';
    b_vals_lpf.name    = 'b-values';
    b_vals_lpf.help    = {'Provide a 1 x N  - matrix with b-values, b-values should appear in the same order as the low- and high-diffusion weighted images were entered. b-values is expected in units of s/mm^2.'
        'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'};
    b_vals_lpf.strtype = 'e';
    b_vals_lpf.num     = [1 Inf];
    b_vals_lpf.val     = {[5 1000 1000 2000]};

    % diffusion directions
    diff_dirs_lpf         = cfg_entry;
    diff_dirs_lpf.tag     = 'diff_dirs_lpf';
    diff_dirs_lpf.name    = 'Diffusion directions';
    diff_dirs_lpf.help    = {'Provide a 3 x N  - matrix with b-vectors, b-vectors should appear in the same order as the low- and high-diffusion weighted images were entered. The b-vectors are dimensionless.'
        'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'};
    diff_dirs_lpf.strtype = 'e';
    diff_dirs_lpf.num     = [3 Inf];
    diff_dirs_lpf.val     = {[1 0 0; 0 1 0; 0 0 1; 0 1/sqrt(2) 1/sqrt(2)]'};

    % in_vols phantom data
    in_vols_lpf_ph         = cfg_files;
    in_vols_lpf_ph.tag     = 'in_vols_lpf_ph';
    in_vols_lpf_ph.name    = 'Phantom DTI images';
    in_vols_lpf_ph.help    = {'Select high- and low-b-value images of the water phantom.'};
    in_vols_lpf_ph.filter  = 'image';
    in_vols_lpf_ph.ufilter = '.*';
    in_vols_lpf_ph.num     = [0 Inf];

    % in_vols diffusion tensor
    in_vols_lpf         = cfg_files;
    in_vols_lpf.tag     = 'in_vols_lpf';
    in_vols_lpf.name    = 'DTI images';
    in_vols_lpf.help    = {'Select high- and low-b-value images of the diffusion data set that has to be correction for LPFs.'};
    in_vols_lpf.filter  = 'image';
    in_vols_lpf.ufilter = '.*';
    in_vols_lpf.num     = [0 Inf];
    
%     diff_lpf         = cfg_exbranch;
%     diff_lpf.tag     = 'diff_lpf';
%     diff_lpf.name    = 'LPF correction of Diffusion Tensor';
%     diff_lpf.val     = {in_vols_lpf in_vols_lpf_ph diff_dirs_lpf b_vals_lpf Nl_lpf perc_lpf epsilon_lpf};
%     diff_lpf.help    = {'Least square tensor estimation using LPF Ellipsoid to correct for perturbation due to gradient inhomogeneities, miscalibration, etc.'
%         'The LPF correction has been first suggested by Bammer et al., MRM 2003.'
%         'The LPF estimation from a water phantom is described in Mohammadi et al., Neuroimage 2012.'
%         'Note that this code is not optimised for RAM space. If the DTI data set is too big, you might get RAM problems.'
%         'Please cite Mohammadi et al., Neuroimage 2012 (doi: 10.1016/j.neuroimage.2011.12.009) when using this code.'};
%     diff_lpf.prog = @local_diff_lpf;
%     diff_lpf.vout = @vout_diff_lpf;
%
%     function out = local_diff_lpf(job)
%     LPFestim_sphharm_PROD(char(job.in_vols_lpf),char(job.in_vols_lpf_ph),job.diff_dirs_lpf,job.b_vals_lpf,job.Nl_lpf, job.perc_lpf, job.epsilon_lpf);
%     if(exist('spm_file','var'))
%         out.FAfiles = spm_file(job.in_vols_lpf{1},'prefix','lpf_'); % to be adjusted
%     else
%         out.FAfiles = my_spm_file(job.in_vols_lpf{1},'prefix','lpf_'); % to be adjusted
%     end
%
%     function dep = vout_diff_lpf(job) %to be adjusted
%     dep(1)            = cfg_dep;
%     dep(1).sname      = 'Diffusion Tensor measures';
%     dep(1).src_output = substruct('.','FAfiles'); % to be adjusted
%     dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); % to be adjusted
    
end

%% Multiple used params
[~, mask, bvals_bvecs, in_vols, ~, p_out, bvals] = tbx_cfg_acid_params;

%% Startup module
startup = tbx_cfg_acid_startup(p_out,bvals_bvecs);

%% PRE-PROCESSING: ECMOCO

% Exbranch 1: Eddy-current and motion correction
ecmoco = tbx_cfg_acid_ecmoco(bvals_bvecs);

% Exbranch 2: Write eddy current and motion correction
ecmoco_write = tbx_cfg_acid_ecmoco_apply;

% Branch: ECMOCO
ecmoco_choice      = cfg_choice;
ecmoco_choice.tag  = 'ecmoco_choice';
ecmoco_choice.name = 'ECMOCO';
ecmoco_choice.help = {
    'Estimate and write ECMOCO'
    'Write ECMOCO'
    }';
ecmoco_choice.values = {ecmoco ecmoco_write};

%% PRE-PROCESSING: msPOAS
poas = tbx_cfg_acid_poas(in_vols, mask, bvals_bvecs);

%% PRE-PROCESSING: HySCO
if dummy_hysco
    hysco_choice = tbx_cfg_acid_hysco(in_vols, dummy_hysco_apply_multi);
end

%% PRE-PROCESSING: Rician bias correction
rbc = tbx_cfg_acid_rbc(in_vols);

rbc_simulation = tbx_cfg_acid_rbc_simulation;

%% Preprocessing
prepro_choice        = cfg_choice;
prepro_choice.tag    = 'prepro_choice';
prepro_choice.name   = 'Pre-processing';
prepro_choice.help   = {
    'Pre-processing of DTI images:'
    '- Displays and creates movies'
    '- Write and read 4d nifti'
    '- Interpolate whole DTI dataset to voxel size of choice'
    '- Combine repetitions'
    '- Rician bias correction'
    '- Eddy current (EC) and Motion correction'
    '- Hyperelastic Susceptibility Artifact Correction (HySCo)'
    '- Position orientation adaptive smoothing (msPOAS)'
    }';
prepro_choice.values = {ecmoco_choice, poas, rbc, hysco_choice};

%% DTI
dummy_dki = 0;
[dti_select,~]  = tbx_cfg_acid_dti(dummy_dki, mask, bvals_bvecs, in_vols);
dti_select.name = 'Diffusion Tensor Imaging (DTI)';

%% Coviper
% Dthr.val = acid_get_defaults('diffusion.dthr');
% diff_coviper = tbx_cfg_acid_coviper(dummy_Tfreiburg, RMatrix, Dthr);

%% UTILITIES: Combine repetitions
% combima = tbx_cfg_acid_combimages;

%% DKI
dummy_dki = 1;
[~,dki_select]  = tbx_cfg_acid_dti(dummy_dki, mask, bvals_bvecs, in_vols);
dki_select.name = 'Diffusion Kurtosis Imaging (DKI)';

%% Diffusion Tensor/Kurtosis imaging
fit_choice        = cfg_choice;
fit_choice.tag    = 'fit_choice';
fit_choice.name   = 'Diffusion tensor/kurtosis imaging';
fit_choice.help   = {
    'Choose the tensor fitting method.'
    }';
fit_choice.values = {dti_select, dki_select};

%% NODDI-DTI
if dummy_noddidti
    acid_noddidti = tbx_cfg_acid_noddidti;
end

%% WMTI-WATSON
acid_wmti_watson = tbx_cfg_acid_wmti_watson;

%% biophysical models
if dummy_noddidti
    biophys_choice        = cfg_choice;
    biophys_choice.tag    = 'biophys_choice';
    biophys_choice.name   = 'Biophysical models';
    biophys_choice.help   = {
        'Choose the biophysical model.'
        }';
    fit_option            = 4;
    diff_BP               = tbx_cfg_acid_gauss_newton_biophysical(fit_option);
    diff_BP.name          = 'NLLS, fit biophysical model directly';
    biophys_choice.values = {acid_noddidti, acid_wmti_watson};
end

%% UTILITIES: Cropping

% Exbranch 1: Create new cropping
crop_new = tbx_cfg_acid_crop;

% Exbranch 2: Apply cropping
crop_apply = tbx_cfg_acid_crop_apply;

% Branch: image cropping
crop      = cfg_choice;
crop.tag  = 'crop';
crop.name = 'Cropping';
crop.help = {
    'New cropping'
    'Apply previous cropping'
    }';
crop.values  = {crop_new, crop_apply};

%% UTILITIES: Resampling
resample = tbx_cfg_acid_resample;

%% UTILITIES: Slice-wise realignment (manual)

% Exbranch 1: manual registration
realign_est = tbx_cfg_acid_realign;

% Exbranch 2: apply manual registration
realign_apply = tbx_cfg_acid_realign_apply;

% Branch: Slice-wise realignment
realign_choice        = cfg_choice;
realign_choice.tag    = 'realign_choice';
realign_choice.name   = 'Slice-wise realignment';
realign_choice.help   = {
    'Slice-wise realignment'
    'Write realigned files'
    }';
realign_choice.values = {realign_est, realign_apply};

%% UTILITIES: Create brain mask
make_brainMSK = tbx_cfg_acid_brainmask;

%% UTILITIES: Reliability masking

% Exbranch 1: determine threshold for reliability masking
relmask_thr = tbx_cfg_acid_relmask_thr;

% Exbranch 2: create reliability masks
relmask_mask = tbx_cfg_acid_relmask_mask;

% Branch: Reliability masking
relmask_choice      = cfg_choice;
relmask_choice.tag  = 'relmask_choice';
relmask_choice.name = 'Reliability masking';
relmask_choice.help = {
    'Determine threshold for masking'
    'Create reliability masks'
    }';
relmask_choice.values = {relmask_thr, relmask_mask};

%% UTILITIES: DWI series browser
browser = tbx_cfg_acid_dwi_series_browser(bvals_bvecs);

%% UTILITIES: DWI series movie
movie = tbx_cfg_acid_dwi_series_movie;

%% UTILITIES: Noise estimation
sigma = tbx_cfg_acid_noise_estimate(bvals);

%% UTILITIES: ROI analysis
roianalysis = tbx_cfg_acid_roi_analysis;

%% UTILITIES: FUSION
fusion = tbx_cfg_acid_fusion;

%% ACID Utilities
image_manipulation_choice = cfg_choice;
image_manipulation_choice.tag    = 'image_manipulation_choice';
image_manipulation_choice.name   = 'Image editing';
image_manipulation_choice.help   = {'ACID Utilities tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
image_manipulation_choice.values = {crop resample, realign_choice, fusion};

qa_choice = cfg_choice;
qa_choice.tag    = 'qa_choice';
qa_choice.name   = 'Quality assessment';
qa_choice.help   = {'ACID Utilities tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
qa_choice.values = {browser, movie};

mask_choice = cfg_choice;
mask_choice.tag    = 'mask_choice';
mask_choice.name   = 'Mask generation';
mask_choice.help   = {'ACID Utilities tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
mask_choice.values = {make_brainMSK, relmask_choice};

utils_choice        = cfg_choice;
utils_choice.tag    = 'utils_choice';
utils_choice.name   = 'Utilities';
utils_choice.help   = {'ACID Utilities tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
utils_choice.values = {image_manipulation_choice, mask_choice, qa_choice, sigma, rbc_simulation, roianalysis};


eddy_call = acid_eddy_call(bvals_bvecs);
wmti_call = acid_wmti_call();

koaya_call = tbx_cfg_acid_rbc_koay(in_vols);

external_preproc_choice        = cfg_choice;
external_preproc_choice.tag    = 'external_preproc_choice';
external_preproc_choice.name   = 'Pre-processing';
external_preproc_choice.help   = {'ACID external tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
external_preproc_choice.values = {eddy_call, koaya_call};

external_biophys_choice        = cfg_choice;
external_biophys_choice.tag    = 'external_biophys_choice';
external_biophys_choice.name   = 'Biophysical Modelling';
external_biophys_choice.help   = {'ACID external tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
external_biophys_choice.values = {wmti_call};

% external_rbc_choice        = cfg_choice;
% external_rbc_choice.tag    = 'external_rbc_choice';
% external_rbc_choice.name   = 'RBC';
% external_rbc_choice.help   = {'ACID external tools. It includes:'
%     '- Displays and creates movies'
%     '- Write and read 4d nifti''Miscellaneous'};
% external_rbc_choice.values = {koaya_call};



external_choice        = cfg_choice;
external_choice.tag    = 'external_choice';
external_choice.name   = 'External tools';
external_choice.help   = {'ACID external tools. It includes:'
    '- Displays and creates movies'
    '- Write and read 4d nifti''Miscellaneous'};
external_choice.values = {external_preproc_choice, external_biophys_choice};


%% DTI artefact correction
acid      = cfg_choice;
acid.tag  = 'dti';
acid.name = 'ACID Toolbox';
acid.help = {
    'Artefact Correction in Diffusion MRI'
    }';

%% ACID batch
if dummy_noddidti
    acid.values = {startup, prepro_choice, fit_choice, biophys_choice, utils_choice, external_choice};
else
    acid.values = {startup, prepro_choice, fit_choice, utils_choice};
end




disp( 'ACID is initialised!                                                                                                ');
fprintf('\n');

end