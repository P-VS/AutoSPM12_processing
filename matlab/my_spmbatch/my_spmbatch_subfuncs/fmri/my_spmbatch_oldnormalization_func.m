function [wfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_oldnormalization_func(ne,funcfile,reffuncfile,ppparams,params,delfiles,keepfiles)

Vfunc = spm_vol(funcfile);
Rfunc= spm_vol(reffuncfile);

for i=1:numel(Vfunc)
    wfuncfiles{i,1} = [funcfile ',' num2str(i)];
end

%% Normalization of the functional scan
if ne==ppparams.echoes(1) || ~isfield(ppparams,'deffile') %Based op OldNorm
    funcnormest.subj.source = {wfuncfiles{1,1}};
    funcnormest.subj.wtsrc = '';
    funcnormest.eoptions.template = {fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii')};
    funcnormest.eoptions.weight = '';
    funcnormest.eoptions.smosrc = 8;
    funcnormest.eoptions.smoref = 0;
    funcnormest.est.eoptions.regtype = 'mni';
    funcnormest.eoptions.cutoff = 25;
    funcnormest.eoptions.nits = 16;
    funcnormest.eoptions.reg = 1;

    spm_run_normalise(funcnormest);

    [fpth,fnm,~] = fileparts(funcfile);

    ppparams.deffile = fullfile(fpth,[fnm '_sn.mat']);
    delfiles{numel(delfiles)+1} = {ppparams.deffile};
end

%% Normalise func

%Write the spatially normalised data

woptions.bb = [-78 -112 -70;78 76 85];
woptions.vox = params.func.normvox;

dt = Vfunc(1).dt;
if dt(1)==spm_type('uint16')
    woptions.interp = 4;
else
    woptions.interp = 1;
end
woptions.prefix = 'w';

funcnormwrite.subj.matname(1) = {ppparams.deffile};
funcnormwrite.subj.resample = wfuncfiles(:,1);
funcnormwrite.roptions.preserve = 0;
funcnormwrite.roptions.bb = woptions.bb;
funcnormwrite.roptions.vox = woptions.vox;
funcnormwrite.roptions.interp = 1;
funcnormwrite.roptions.wrap = [0 0 0];
funcnormwrite.roptions.prefix = woptions.prefix;

spm_run_normalise(funcnormwrite);

ppparams.func(ne).wfuncfile = ['w' ppparams.func(ne).prefix ppparams.func(ne).funcfile];
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ne).wfuncfile)};

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(ne).wfuncfile));
wfuncdat = spm_read_vols(Vfunc);