function [ppparams,keepfiles,delfiles] = fmri_ica_aroma(ppparams,keepfiles,delfiles)
%FMRI_ICA_AROMA Performs ICA-AROMA

%% ICA step

% Initialize parameters and run ICA:

% Get t_r
jsondat = jsondecode(fileread(ppparams.funcjsonfile));
t_r = jsondat.("RepetitionTime");

curdir = pwd;

do_ica(ppparams.fmask, t_r, ppparams);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'input_spatial_ica.m')};
delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'ica_dir')};

cd(curdir)

%% AROMA (IC classification)

ica_dir = fullfile(ppparams.ppfuncdir, 'ica_dir'); % path to ICs computed at previous step

noiseICdata = ica_aroma_classification(ppparams, ppparams.fmask, ica_dir, t_r);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'headdata.nii')};

%% Denoising 

ppparams.rp_file = spm_file(ppparams.rp_file, 'prefix','ica-aroma_','ext','.txt');

writematrix(noiseICdata,ppparams.rp_file,'Delimiter','tab');

end

