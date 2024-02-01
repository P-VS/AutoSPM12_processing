function AutoSPMprocessing_fmri

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

params.spm_path = '/Users/accurad/Library/CloudStorage/OneDrive-Personal/Matlab/spm12';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/rTMS-fMRI_Interleaved/Data';

first_sub = 1;
last_sub = 1;
sublist = [3]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0
 
params.task = {'semantic'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

params.analysisname = 'test';

params.use_parallel = true; 
params.maxprocesses = 4; %Best not too high to avoid memory problems
params.keeplogs = false;

%% fMRI data parameters
    params.preprocfmridir = 'preproc_func_test'; %directory with the preprocessed fMRI data
    params.fmri_prefix = 'sware'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']
    params.fmri_endfix = 'bold';
    
    params.dummytime = 0;
    
    params.multi_echo = false;
    params.echoes = [1,2,3]; %list of echoes for ME-fMRI used as sessions in SPM
    params.use_echoes_as_sessions = false; %use each echo series as a session in SPM

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

    params.contrast(1).conditions = {'semantic'};
    params.contrast(1).vector = [1];
    
    params.contrast(2).conditions = {'semantic'};
    params.contrast(2).vector = [-1];
    
    %params.contrast(1).conditions = {'episodic','semantic'};
    %params.contrast(1).vector = [1,-1];
    
    %params.contrast(2).conditions = {'episodic','semantic'};
    %params.contrast(2).vector = [-1,1];
    
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

fprintf('\nDone with processing the data\n')

end