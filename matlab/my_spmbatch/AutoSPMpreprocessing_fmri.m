function AutoSPMpreprocessing_fmri

%Script to do the auto preprocessing in SPM12
%
%Preparation:
% Convert the DICOM files into nifti using dcm2niix in MROCroGL
% For the anatomical scans, set 'Crop 3D Images' on
%
%* Organise the data in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1)
%                   Files: sub-##_T1w.nii and sub-##_T1w.json
%                -func: containes the fmri data
%                   Files: sub-##_task-..._bold.nii and sub-##_task-..._bold.json
%                -fmap: containnes the gradient pololarity (blip-up/down) filpt data or the fieldmap scans
%                   Files in case of inverted gradient polarity: sub-##_dir-pi_epi.nii and sub-##_dir-pi_epi.json
%                   Files in case of fieldmap scans: (image 1 in file is amplitude, image 2 in file is phase)
%                          sub-##_fmap_echo-1.nii and sub-##_fmap_echo-1.json
%                          sub-##_fmap_echo-2.nii and sub-##_fmap_echo-2.json
    
%* IMPORTANT: !! Look at your data before starting any (pre)processing. Losing time in trying to process bad data makes no sense !!

%Script written by dr. Peter Van Schuerbeek (Radiology UZ Brussel)

%% Give path to SPM12 and GroupICA

params.spm_path = '/Users/accurad/Library/CloudStorage/OneDrive-Personal/Matlab/spm12';
params.GroupICAT_path = '/Users/accurad/Library/CloudStorage/OneDrive-Personal/Matlab/GroupICATv40c';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/rTMS-fMRI_Interleaved/Data';

sublist = [3];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0

params.save_folder = 'preproc_func_test';

task ={'semantic'};

params.use_parallel = false; 
params.maxprocesses = 4; %Best not too high to avoid memory problems
params.keeplogs = false;

params.save_intermediate_results = false; 

params.reorient = true; % align data with MNI template to improve normalization and segmentation

%% Preprocessing anatomical data

    params.preprocess_anatomical = true;

    % Normalization
    params.anat.do_normalization = true;
    params.anat.normvox = [1.5 1.5 1.5];

    %% When not doing VBM -> Use default in SPM
        % Segmentation
        params.anat.do_segmentation = true;
    
%% Preprocessing functional data

    params.preprocess_functional = false;

    % Remove the dummy scans n_dummy_scans = floor(dummytime/TR)
    params.func.dummytime = 6; %time in seconds

    % Realignnment (motion correction)
    params.func.do_realignment = true;
    
    % Geometric correction
    % fieldmap or pepolar Only 1 can be true, the other should be false.
    params.func.fieldmap = false;
    params.func.pepolar = false;

    % Slice time correction
    params.func.do_slicetime = true;
        
    % If ME-fMRI, combine the multiple eccho images
    params.func.meepi = false;
    params.func.echoes = [1]; %the numbers of the echoes in ME-fMRI. 
    params.func.combination = 'none'; 
    %none: all echoes are preprocessed separatly
    %average: The combination is the average of the multiple echo images
    %TE_weighted: The combination is done wi=TEi or 
    %T2_weighted: dynamic T2* weighted combination
    %dyn_T2star: dynamic T2* mapping
    %see Heunis et al. 2021. The effects of multi-echo fMRI combination and rapid T2*-mapping on offline and real-time BOLD sensitivity. NeuroImage 238, 118244

    % Normalization
    params.func.do_normalization = true;
    params.func.normvox = [1.5 1.5 1.5];

    % Smoothing
    params.func.do_smoothing = true;
    params.func.smoothfwhm = 6;

%% Denoising

    params.do_denoising = false; 

        %if do_denoising = true and preprocess_functional = false -> only denoising, other preprocessing already done
        params.denoise.prefix = 'swaure';

    params.denoise.meepi = false;
    params.denoise.echoes = [1]; %the numbers of the echoes in ME-fMRI. 
    params.denoise.mecombined = false; %if ME-fMRI, where the echoes combined?
    
    % Band-pass filtering
    params.denoise.do_bpfilter = false;
    params.denoise.bpfilter = [0.008 Inf]; %no highpass filter is first 0, no lowpass filter is last Inf, default is [0.008 0.1]

    % Extend motion regressors with derivatives and squared regressors
    params.denoise.do_mot_derivatives = true; %derivatives+squares (24 regressors)

    % aCompCor
    params.denoise.do_aCompCor = true;
    params.denoise.Ncomponents = 5; %if in range [0 1] then the number of aCompCor components is equal to the number of components that explain the specified percentage of variation in the signal

    % ICA-AROMA
    params.denoise.do_ICA_AROMA = true;

    % Noise regression / remove ICA-AROMA noise components
    params.denoise.do_noiseregression = true;

    % Prepare data for DENN denoising in python
    params.denoise.do_DENN = false;
    
%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------

restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if exist(params.GroupICAT_path,'dir'), addpath(genpath(params.GroupICAT_path)); end
if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end

fprintf('Start with preprocessing \n')

curdir = pwd;

warnstate = warning;
warning off;

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

spm('defaults', 'FMRI');

my_spmbatch_start_fmripreprocessing(sublist,nsessions,task,datpath,params)

spm_figure('close',allchild(0));

cd(curdir)

fprintf('\nDone\n')

end