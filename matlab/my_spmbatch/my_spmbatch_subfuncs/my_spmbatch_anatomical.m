function [delfiles,keepfiles] = my_spmbatch_anatomical(sub,ses,datpath,params)

ppparams.sub = sub;
ppparams.ses = ses;

ppparams.substring = ['sub-' num2str(sub,'%02d')];

if ~isfolder(fullfile(datpath,ppparams.substring))
    ppparams.substring = ['sub-' num2str(sub,'%03d')];
end

ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%03d')]);

if ~isfolder(ppparams.subpath)
    ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%02d')]);
end

ppparams.subanatdir = fullfile(ppparams.subpath,'anat');
ppparams.preproc_anat = fullfile(ppparams.subpath,'preproc_anat');

delfiles = {};
keepfiles = {};

[delfiles,keepfiles] = my_spmbatch_preprocess_anat(ppparams,params,delfiles,keepfiles);