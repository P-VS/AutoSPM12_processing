function [ppparams,keepfiles] = my_spmbatch_acompcor(funcdat,ppparams,params,keepfiles)

if isfield(ppparams,'der_file')
    confounds = load(ppparams.der_file);
elseif isfield(ppparams,'rp_file')
    confounds = load(ppparams.rp_file);
else
    confounds = [];
end

GM = spm_vol(ppparams.fc1im);
WM = spm_vol(ppparams.fc2im);
CSF = spm_vol(ppparams.fc3im);

gmdat = spm_read_vols(GM);
wmdat = spm_read_vols(WM);
csfdat = spm_read_vols(CSF);

braindat = gmdat+wmdat;
braindat(braindat<0.2)=0;
braindat(braindat>0.0)=1;

csfdat(braindat>0.0)=0;
csfdat(csfdat<0.2)=0;
csfdat(csfdat>0.0)=1;

acc_confounds = fmri_acompcor(funcdat(:,:,:,:),{csfdat},params.denoise.Ncomponents,'confounds',confounds,'filter',[],'PolOrder',1);

accf = ['acc_' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile];
accname = split(accf,'_echo-');
accname = split(accf,'.nii');
ppparams.acc_file = fullfile(ppparams.subfuncdir,[accname{1} '.txt']);

writematrix(acc_confounds,ppparams.acc_file,'Delimiter','tab');

keepfiles{numel(keepfiles)+1} = {ppparams.acc_file};

clear csfdat braindat gmdat wmdat