function [wfuncdat,ppparams,keepfiles] = my_spmbatch_noiseregression(wfuncdat,ne,ppparams,params,keepfiles)

fprintf('Start noise regression \n')

denprefix = 'd';

if ~isfield(ppparams,'noiseregresssor')
    if params.denoise.do_ICA_AROMA && isfield(ppparams,'ica_file')
        confounds = load(ppparams.ica_file);
    else
        if params.denoise.do_mot_derivatives && isfield(ppparams,'der_file')
            confounds = load(ppparams.der_file);
        elseif isfield(ppparams,'rp_file')
            confounds = load(ppparams.rp_file);
        else
            confounds = [];

            denprefix = 'f';
        end
    
        if params.denoise.do_aCompCor && isfield(ppparams,'acc_file')
            acc_confounds = load(ppparams.acc_file);

            confounds = cat(2,confounds,acc_confounds);
        end
    end
else
    confounds = ppparams.noiseregresssor;
end

if params.denoise.do_bpfilter
    jsondat = fileread(ppparams.funcjsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    bpfilter = [tr ppparams.bpfilter(1:2)];
else
    bpfilter = [];
end

if exist('wfuncdat','var')
    s = size(wfuncdat);
    wfuncdat = reshape(wfuncdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

    [wfuncdat,~] = fmri_cleaning(wfuncdat(:,:),1,bpfilter,confounds,[],'restoremean','on');
else
    [wfuncdat,~] = fmri_cleaning(fullfile(ppparams.ppfuncdir,ppparams.func(ne).sfuncfile),1,bpfilter,confounds,[],'restoremean','on');
end

wfuncdat = reshape(wfuncdat(:,:),s);

Vfunc = spm_vol(fullfile(ppparams.ppfuncdir,ppparams.func(ne).sfuncfile));

for k=1:numel(Vfunc)
    Vfunc(k).fname = fullfile(ppparams.ppfuncdir,[denprefix ppparams.func(ne).sfuncfile]);
    Vfunc(k).descrip = 'my_spmbatch - denoise';
    Vfunc(k).pinfo = [1,0,0];
    Vfunc(k).n = [k 1];
end

Vfunc = myspm_write_vol_4d(Vfunc,wfuncdat);

ppparams.func(ne).dfuncfile = [denprefix ppparams.func(ne).sfuncfile];
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.ppfuncdir,[denprefix ppparams.func(ne).sfuncfile])};

fprintf('Done noise regression \n')