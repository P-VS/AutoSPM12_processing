function AutoSPM1stlevel_fmri

%Script to do the auto 1st level fMRI processing in SPM12
%
%Preparation:
%* Organise the data in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1)
%                -func: containes the fmri data
%            
%* Make sure all data is preprocessed
%    - preproc_anat
%    - preproc_func
    
%* IMPORTANT: !! Look at your preprocessed data before starting any analysis. It makes no sense to lose time in trying to process bad data !!

%% Give path to SPM12

clear 'all'

params.spm_path = '/Users/accurad/Library/CloudStorage/OneDrive-Personal/Matlab/spm12';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/asl_bold/ASL_fingertapping';

first_sub = 1;
last_sub = 1;
sublist = [1:13]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0
 
params.task = {'bilateralfingertapping'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

params.analysisname = 'func_meica';
params.modality = 'fmri'; %'fmri' or 'fasl'

params.use_parallel = false; 
params.maxprocesses = 2; %Best not too high to avoid memory problems
params.loadmaxvols = 100; %to reduce memory load, the preprocessing can be split in smaller blocks (default = 100)
params.keeplogs = false;

%% fMRI data parameters
    params.preprocfmridir = 'preproc_func_meica'; %directory with the preprocessed fMRI data
    params.fmri_prefix = 'swacdfre'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']
    
    params.dummytime = 0; %only if the timings in the _events.tsv file should be corrected for dummy scans
    
    %In case of multiple runs in the same session exist
    params.func.mruns = false; %true if run number is in filename
    params.func.runs = [1]; %the index of the runs (in filenames run-(index))
    params.func.use_runs = 'together'; % 'separately' or 'together' (this parameter is ignored if mruns is false)
    %'separately': a separate analysis is done per run
    %'together': all runs are combined in 1 analysis
    
    % For ME-fMRI
    params.func.meepi = true; %true if echo number is in filename
    params.func.echoes = [1,2,3,4]; %the index of echoes in ME-fMRI used in the analysis. If meepi=false, echoes=[1]. 

%% SPM first level analysis parameters
    params.confounds_prefix = 'rp_e'; %confounds file of form [confounds_prefix 'sub-ii_task-... .txt']
    params.add_regressors = false;
    params.use_ownmask = false;
    params.model_serial_correlations = 'AR(1)';
    params.hpf = 128; %default 128
    
    %contrast(i) is structure with fields
    %   conditions: conditions to compare
    %   vector: contrast weight vector
    %e.g A contrast between conditions is given as
    %   contrast(i).conditions={'condition 1','condition 2'};
    %   contrast(i).vector=[1 -1];

    params.contrast(1).conditions = {'Finger'};
    params.contrast(1).vector = [1];

    params.contrast(2).conditions = {'Finger'};
    params.contrast(2).vector = [-1];

    %params.contrast(1).conditions = {'episodic','semantic'};
    %params.contrast(1).vector = [1,-1];
    
    %params.contrast(2).conditions = {'episodic','semantic'};
    %params.contrast(2).vector = [-1,1];
    
    %params.contrast(3).conditions = {'episodic','semantic'};
    %params.contrast(3).vector = [1,1];
    
    %params.contrast(4).conditions = {'episodic','semantic'};
    %params.contrast(4).vector = [-1,-1];
    
    %params.contrast(1).conditions = {'sad','neutral'};
    %params.contrast(1).vector = [1,-1];
    
    %params.contrast(2).conditions = {'sad','neutral'};
    %params.contrast(2).vector = [-1,1];
    
    %params.contrast(3).conditions = {'happy','neutral'};
    %params.contrast(3).vector = [1,-1];
    
    %params.contrast(4).conditions = {'happy','neutral'};
    %params.contrast(4).vector = [-1,1];
    
    %params.contrast(5).conditions = {'sad','happy'};
    %params.contrast(5).vector = [1,-1];
    
    %params.contrast(6).conditions = {'sad','happy'};
    %params.contrast(6).vector = [-1,1];

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%--------------------------------------------------------------------------------

restoredefaultpath

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end

old_spm_read_vols_file=fullfile(spm('Dir'),'spm_read_vols.m');
new_spm_read_vols_file=fullfile(spm('Dir'),'old_spm_read_vols.m');

if isfile(old_spm_read_vols_file), movefile(old_spm_read_vols_file,new_spm_read_vols_file); end
  
fprintf('Start with processing the data\n')

warnstate = warning;
warning off;

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

spm('defaults', 'FMRI');

curdir = pwd;

my_spmbatch_start_fmriprocessing(sublist,nsessions,datpath,params);

spm_figure('close',allchild(0));

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

fprintf('\nDone with processing the data\n')

clear 'all'

end