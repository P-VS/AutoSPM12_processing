function [delfiles,keepfiles] = my_spmbatch_functional(sub,ses,task,datpath,params)

ppparams.testscript = false;

if ~params.func.meepi
    params.func.echoes = [1];
end

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.task = task;
ppparams.use_parallel = params.use_parallel;

ppparams.substring = ['sub-' num2str(sub,'%02d')];

if ~isfolder(fullfile(datpath,ppparams.substring))
    ppparams.substring = ['sub-' num2str(sub,'%03d')];
end

ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%03d')]);

if ~isfolder(ppparams.subpath)
    ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%02d')]);
end

ppparams.subfmridir = fullfile(ppparams.subpath,'func');

delfiles = {};
keepfiles = {};

[delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles);