function [wfuncdat,funcfile,keepfiles] = my_spmbatch_noiseregression(subfmridir,substring,task,funcfile,rp_file,wfuncdat,keepfiles)

fprintf('Start noise regression \n')

if exist(rp_file)
    confounds = load(rp_file);
else
    confounds = [];
end

if params.do_bpfilter
    if params.nechoes==1
        funcjsonfile = fullfile(subfmridir,[substring '_task-' task '_bold.json']);
    else
        funcjsonfile = fullfile(subfmridir,[substring '_task-' task '_bold_e1.json']);
    end

    jsondat = fileread(funcjsonfile);
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
    [wfuncdat,~] = fmri_cleaning(funcfile,1,bpfilter,confounds,[],'restoremean','on');
end

wfuncdat = reshape(wfuncdat(:,:),s);

Vfunc = spm_vol(funcfile);

for k=1:numel(Vfunc)
    Vfunc(k).fname = spm_file(funcfile, 'prefix','d');
    Vfunc(k).descrip = 'my_spmbatch - denoise';
    Vfunc(k).n = [k 1];
    Vfunc(k) = spm_create_vol(Vfunc(k));
    Vfunc(k) = spm_write_vol(Vfunc(k),wfuncdat(:,:,:,k));
end

funcfile = spm_file(funcfile, 'prefix','d');
keepfiles{numel(keepfiles)+1} = {funcfile};

fprintf('Done noise regression \n')