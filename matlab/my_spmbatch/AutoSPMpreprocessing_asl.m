function AutoSPMpreprocessing_asl

%Script to do the auto ASL preprocessing in SPM12
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
%                -perf: containes the ASL data
%                   Files: sub-##_asl_m0scan.nii and sub-##_asl_m0scan.json
%                          sub-##_asl_deltam.nii and sub-##_asl_deltam.json
%                          sub-##_asl_cbf.nii and sub-##_asl_cbf.json
%    
%* IMPORTANT: !! Look at your data before starting any (pre)processing. Losing time in trying to process bad data makes no sense !!

%Script written by dr. Peter Van Schuerbeek (Radiology UZ Brussel)

%% Give path to SPM12

params.spm_path = '/Users/accurad/Library/CloudStorage/OneDrive-Personal/Matlab/spm12';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_studie_Linde/Data';

sublist = [5,8:11,13:14,16:20];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [2]; %nsessions>0

params.multiple_scans = false;
params.num_scans = [1]; %File names are sub##_asl-#_...

params.save_folder = 'preproc_aslge_gm';

params.use_parallel = true; 
params.maxprocesses = 4; %Best not too high to avoid memory problems
params.keeplogs = false;

params.save_intermediate_results = false; 

params.reorient = true; % align data with MNI template to improve normalization

%% In case data already processsed at the GE scanner (PCASL scan)
    
    params.asl.preprocessed_ge = true;
    params.aslge.cbfmap_present = false;

%% In case data already processsed at the GE scanner (PCASL scan)

    params.aslge.do_cbfmapping = true;
    params.aslge.T1correctionM0 = 'average_GM'; %T1 correction of M0scan: 'tisssue_maps', 'T1_map', 'average_GM', 'average_WM'

%% Segmentation

    params.asl.do_segmentation = false;

%% Normalization

    params.asl.do_normalization = true;
    params.asl.normvox = [2 2 2];
  
%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------

restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

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

my_spmbatch_start_aslpreprocessing(sublist,nsessions,datpath,params)

spm_figure('close',allchild(0));

cd(curdir)

fprintf('\nDone\n')

end