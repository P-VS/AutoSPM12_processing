function AutoSPMpreprocessing

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

warnstate = warning;
warning off;

%spm_figure('Clear','Interactive')

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

%%
%%Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data';

sublist = [1,2,4:9];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0

params.save_folder = 'preproc_func_dyn-t2star';

task ={'ME-EmoFaces'};

params.meepi = true;
params.echoes = [1,2,3]; %number of echoes for ME-fMRI. 
params.combination = 'dyn_T2star'; 
%none: all echoes are preprocessed separatly
%average: The combination is the average of the multiple echo images
%TE_weighted: The combination is done wi=TEi or 
%T2_fit: dynamic T2* weighted combination
%dyn_T2star: dynamic T2* mapping
%see Heunis et al. 2021. The effects of multi-echo fMRI combination and rapid T2*-mapping on offline and real-time BOLD sensitivity. NeuroImage 238, 118244

params.dummytime = 8; %time in seconds

params.reorient = true;

params.do_segmentation = false;
params.do_slicetime = true;

%fieldmap or pepolar should be true, the other should be false
params.fieldmap = false;
params.pepolar = true;

params.do_realignment = true;

params.do_normalization = true;
params.do_smoothing = true;

params.normvox = [1.5 1.5 1.5];
params.smoothfwhm = 6;

params.do_bpfilter = false;
params.do_mot_derivatives = false; %derivatives+squares (24 regressors)
params.do_aCompCor = false;
params.do_noiseregression = false;
params.do_ICA_AROMA = false;

params.bpfilter = [0.008 Inf]; %no highpass filter is first 0, no lowpass filter is last Inf, default is [0.008 0.1]
params.Ncomponents = 5; %if in range [0 1] then the number of aCompCor components is equal to the number of components that explain the specified percentage of variatiion in the signal

use_parallel = true;
save_intermediate_results = false;

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------
fprintf('Start with preprocessing \n')

curdir = pwd;

spm('defaults', 'FMRI');

if use_parallel
    my_spmbatch_parallel(sublist,nsessions,task,datpath,params,save_intermediate_results)
else
    my_spmbatch_noparallel(sublist,nsessions,task,datpath,params,save_intermediate_results)
end

cd(curdir)

fprintf('\nDone\n')

end