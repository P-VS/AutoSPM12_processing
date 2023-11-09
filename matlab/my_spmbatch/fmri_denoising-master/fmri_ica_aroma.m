function [ppparams,keepfiles,delfiles] = fmri_ica_aroma(tedata,ppparams,keepfiles,delfiles)
%FMRI_ICA_AROMA Performs ICA-AROMA

%% Make masks

% Functional mask
func_mask = my_spmbatch_mask(tedata{ppparams.echoes(1)}.wfuncdat);

Vfunc = tedata{ppparams.echoes(1)}.Vfunc;
Vfuncmask = Vfunc(1);
Vfuncmask.fname = fullfile(ppparams.ppfuncdir, 'funcmask.nii'); % Change name to contain subjectID
Vfuncmask.descrip = 'funcmask';
Vfuncmask = rmfield(Vfuncmask, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
spm_write_vol(Vfuncmask, func_mask);

delfiles{numel(delfiles)+1} = {Vfuncmask.fname};   

%% ICA step

% Initialize parameters and run ICA:

% Get t_r
jsondat = jsondecode(fileread(ppparams.funcjsonfile));
t_r = jsondat.("RepetitionTime");

curdir = pwd;

do_ica(Vfuncmask.fname, tedata, numel(tedata), t_r, ppparams);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'input_spatial_ica.m')};
delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'ica_dir')};

cd(curdir)

%% AROMA (IC classification)

ica_dir = fullfile(ppparams.ppfuncdir, 'ica_dir'); % path to ICs computed at previous step

noiseICdata = ica_aroma_classification(ppparams, Vfuncmask.fname, ica_dir, t_r);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir, 'headdata.nii')};

%% Denoising 

ppparams.rp_file = spm_file(ppparams.rp_file, 'prefix','ica-aroma_','ext','.txt');

writematrix(noiseICdata,ppparams.rp_file,'Delimiter','tab');

end

