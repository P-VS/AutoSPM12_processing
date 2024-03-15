function [ppparams,keepfiles] = my_spmbatch_acompcor(wfuncdat,ppparams,params,keepfiles)

if isfield(ppparams,'der_file')
    confounds = load(ppparams.der_file);
elseif isfield(ppparams,'rp_file')
    confounds = load(ppparams.rp_file);
else
    confounds = [];
end

GM = spm_vol(ppparams.wc1im);
WM = spm_vol(ppparams.wc2im);
CSF = spm_vol(ppparams.wc3im);

gmdat = spm_read_vols(GM);
wmdat = spm_read_vols(WM);
csfdat = spm_read_vols(CSF);

braindat = gmdat+wmdat;
braindat(braindat<0.2)=0;
braindat(braindat>0.0)=1;

csfdat(braindat>0.0)=0;
csfdat(csfdat<0.2)=0;
csfdat(csfdat>0.0)=1;

if params.denoise.do_bpfilter
    jsondat = fileread(ppparams.funcjsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    bpfilter = [tr ppparams.bpfilter(1:2)];
else
    bpfilter = [];
end

acc_confounds = fmri_acompcor(wfuncdat(:,:,:,:),{csfdat},ppparams.Ncomponents,'confounds',confounds,'filter',bpfilter,'PolOrder',1);

ppparams.acc_file = spm_file(fullfile(ppparams.ppfuncdir,ppparams.func(1).wfuncfile), 'prefix','acc_','ext','.txt');

writematrix(acc_confounds,ppparams.acc_file,'Delimiter','tab');

keepfiles{numel(keepfiles)+1} = {ppparams.acc_file};
