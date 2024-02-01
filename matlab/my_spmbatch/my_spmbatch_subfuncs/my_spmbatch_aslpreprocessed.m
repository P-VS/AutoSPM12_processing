function [delfiles,keepfiles] = my_spmbatch_aslpreprocessed(sub,ses,datpath,params)

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

ppparams.subasldir = fullfile(ppparams.subpath,'perf');

delfiles = {};
keepfiles = {};

[delfiles,keepfiles] = my_spmbatch_preprocasl_processed(ppparams,params,delfiles,keepfiles);