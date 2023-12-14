function [ppparams,keepfiles] = my_spmbatch_acompcor(wfuncdat,ppparams,params,keepfiles)

if exist(ppparams.rp_file,'file')
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

if exist(ppparams.rp_file)
    confounds = cat(2,confounds,acc_confounds);

    ppparams.rp_file = spm_file(ppparams.rp_file, 'prefix','acc_','ext','.txt');

    writematrix(confounds,ppparams.rp_file,'Delimiter','tab');
else
    ppparams.rp_file = spm_file(ppparams.funcfile, 'prefix','acc_','ext','.txt');

    writematrix(acc_confounds,ppparams.rp_file,'Delimiter','tab');
end

keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
