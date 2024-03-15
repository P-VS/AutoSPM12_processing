function [ppparams,keepfiles,delfiles] = fmri_ica_aroma(ppparams,keepfiles,delfiles)
%FMRI_ICA_AROMA Performs ICA-AROMA

%% ICA step

% Initialize parameters and run ICA:

% Get t_r
jsondat = jsondecode(fileread(ppparams.funcjsonfile));
t_r = jsondat.("RepetitionTime");

curdir = pwd;


ica_dir = fullfile(ppparams.ppfuncdir, 'ica_dir'); % path to ICs computed at previous step

if ~exist(ica_dir,"dir")
    fprintf('Start ICA \n') 

    do_ica(ppparams.fmask, t_r, ppparams); 

    delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'input_spatial_ica.m')};
end

cd(curdir)

%% AROMA (IC classification)

fprintf('Start ICA classification\n')

noiseICdata = ica_aroma_classification(ppparams, ppparams.fmask, ica_dir, t_r);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'headdata.nii')};

%% Denoising 

ppparams.ica_file = spm_file(fullfile(ppparams.ppfuncdir,ppparams.func(1).sfuncfile), 'prefix','ica_','ext','.txt');

writematrix(noiseICdata,ppparams.ica_file,'Delimiter','tab');

end

