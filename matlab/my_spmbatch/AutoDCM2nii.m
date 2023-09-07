function AutoDCM2nii

%Script to do the auto convertion from DCM file into nii files
%
%* The nii data will be saved in BIDS format
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

params.datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data';

sublist = [3];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]

%Give per sequence
%folder: substructure starting from sub-##
%seqtype: (anat, func, fmap, dti, asl) will be used to organise the data in folders
%name: name of the sequence used in sub-##_... .nii to add series number to the name use %s, to add echo number use %e

params.mridata(1).folder = 'DCM/3DT1';
params.mridata(1).seqtype = 'anat';
params.mridata(1).name = 'T1w';
params.mridata(1).session = 1;

params.mridata(2).folder = 'DCM/ME-fMRI_PI';
params.mridata(2).seqtype = 'fmap';
params.mridata(2).name = 'dir-pi_epi';
params.mridata(2).session = 1;

params.mridata(3).folder = 'DCM/ME-fMRI_EFT';
params.mridata(3).seqtype = 'func';
params.mridata(3).name = 'task-ME-EFT_bold';
params.mridata(3).session = 1;

use_parallel = true; %only possible when parallel toolbox is installed

%Be careful with changing the code below this line!
%--------------------------------------------------------------------------------
fprintf('Start with converting the data\n')

spm('defaults', 'FMRI');

curdir = pwd;

if use_parallel

    my_spmbatch_paralleldcm2nii(sublist,params)

else

    my_spmbatch_noparalleldcm2nii(sublist,params)

end

cd(curdir)

fprintf('\nDone with converting the data\n')

end