function [wfuncdat,ppparams,keepfiles] = my_spmbatch_noiseregression(wfuncdat,ne,ppparams,params,keepfiles)

fprintf('Start noise regression \n')

if exist(ppparams.rp_file)
    confounds = load(ppparams.rp_file);
else
    confounds = [];
end

if params.do_bpfilter
    jsondat = fileread(ppparams.funcjsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    bpfilter = [tr params.bpfilter(1:2)];
else
    bpfilter = [];
end

if exist('wfuncdat','var')
    s = size(wfuncdat);
    wfuncdat = reshape(wfuncdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

    [wfuncdat,~] = fmri_cleaning(wfuncdat(:,:),1,bpfilter,confounds,[],'restoremean','on');
else
    [wfuncdat,~] = fmri_cleaning(ppparams.funcfile,1,bpfilter,confounds,[],'restoremean','on');
end

wfuncdat = reshape(wfuncdat(:,:),s);

Vfunc = spm_vol(ppparams.funcfile{ne});

for k=1:numel(Vfunc)
    Vfunc(k).fname = spm_file(ppparams.funcfile{ne}, 'prefix','d');
    Vfunc(k).descrip = 'my_spmbatch - denoise';
    Vfunc(k).pinfo = [1,0,0];
    Vfunc(k).n = [k 1];
end

Vfunc = myspm_write_vol_4d(Vfunc,wfuncdat);

ppparams.funcfile{ne} = spm_file(ppparams.funcfile{ne}, 'prefix','d');
keepfiles{numel(keepfiles)+1} = {ppparams.funcfile{ne}};

fprintf('Done noise regression \n')