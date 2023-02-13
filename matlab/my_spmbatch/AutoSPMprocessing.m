function AutoSPMprocessing

%Script to do the auto processing in SPM12
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
%    - funcmask.nii
    
%* IMPORTANT: !! Look at your preprocessed data before starting any analysis. It makes no sense to lose time in trying to process bad data !!

%﻿Give the basic input information of your data

params.datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data';
params.analysisname = 'SE_EFT_ICA-AROMA';

first_sub = 1;
last_sub = 1;
sublist = [2:7]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0

task = {'SE-EFT'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

params.preprocfmridir = 'preproc_func_ICA-AROMA';
params.fmri_prefix = 'swuae'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']
params.fmri_endfix = 'bold_regfilt';

params.dummytime = 0;
params.multi_echo=false;

params.confounds_prefix = 'rp_ae'; %confounds file of form [confounds_prefix 'sub-ii_task-... .txt']
params.add_regressors = false;
params.use_ownmask = false;
params.model_serial_correlations = 'AR(1)';
params.hpf = 200; %default 128

%contrast(i) is structure with fields
%   conditions: conditions to compare
%   vector: contrast weight vector
%e.g A contrast between conditions is given as
%   contrast(i).conditions={'condition 1','condition 2'};
%   contrast(i).vector=[1 -1];

params.contrast(1).conditions = {'episodic','semantic'};
params.contrast(1).vector = [1,-1];

params.contrast(2).conditions = {'episodic','semantic'};
params.contrast(2).vector = [-1,1];

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

use_parallel = true;

%Be careful with changing the code below this line!
%--------------------------------------------------------------------------------
fprintf('Start with processing the data\n')

spm('defaults', 'FMRI');

curdir = pwd;

if use_parallel
    datlist = zeros(numel(sublist)*numel(nsessions),2);

    dpos = 1;
    for i = 1:numel(sublist)
        for j = 1:numel(nsessions)
            datlist(dpos,1) = sublist(i);
            datlist(dpos,2) = nsessions(j);

            dpos = dpos+1;
        end
    end

    pa=parpool(min([5,numel(datlist(:,1))])); %25 is the maximum number of workers allowed in the 'local' profile while 10 is set to avoid memory issues on my computer
    parfor i = 1:numel(datlist(:,1))
                
        %% make batch
        matlabbatch = my_spmprocessingbatch(datlist(i,1),datlist(i,2),task,params);
    
        spm_jobman('run', matlabbatch);
    end
    delete(pa)
else
    for i = sublist
        for j = nsessions
            matlabbatch = my_spmprocessingbatch(i,j,task,params);

            spm_jobman('run', matlabbatch);
        end
    end
end

fprintf('\nDone with processing the data\n')
close all

end