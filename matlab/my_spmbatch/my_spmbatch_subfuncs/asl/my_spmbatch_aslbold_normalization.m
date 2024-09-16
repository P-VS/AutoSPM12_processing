function [ppparams,delfiles,keepfiles] = my_spmbatch_aslbold_normalization(ppparams,params,delfiles,keepfiles)

%% Normalization of the functional scan based op OldNorm
aslnormest.subj.source = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile])};
aslnormest.subj.wtsrc = '';
aslnormest.eoptions.template = {fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii')};
aslnormest.eoptions.weight = '';
aslnormest.eoptions.smosrc = 8;
aslnormest.eoptions.smoref = 0;
aslnormest.est.eoptions.regtype = 'mni';
aslnormest.eoptions.cutoff = 25;
aslnormest.eoptions.nits = 16;
aslnormest.eoptions.reg = 1;

spm_run_normalise(aslnormest);

fnm = split(ppparams.perf(1).m0scanfile,'.nii');

ppparams.deffile = fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix fnm{1} '_sn.mat']);
delfiles{numel(delfiles)+1} = {ppparams.deffile};

%% Normalise CBF data

Vcbf = spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile]));

for i=1:numel(Vcbf)
    wcbffiles{i,1} = [Vcbf(i).fname ',' num2str(i)];
end

%Write the spatially normalised data

woptions.bb = [-78 -112 -70;78 76 85];
woptions.vox = params.func.normvox;

dt = Vcbf(1).dt;
if dt(1)==spm_type('uint16')
    woptions.interp = 4;
else
    woptions.interp = 1;
end
woptions.prefix = 'w';

aslnormwrite.subj.matname(1) = {ppparams.deffile};
aslnormwrite.subj.resample = wcbffiles(:,1);
aslnormwrite.roptions.preserve = 0;
aslnormwrite.roptions.bb = woptions.bb;
aslnormwrite.roptions.vox = woptions.vox;
aslnormwrite.roptions.interp = 1;
aslnormwrite.roptions.wrap = [0 0 0];
aslnormwrite.roptions.prefix = woptions.prefix;

spm_run_normalise(aslnormwrite);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['w' ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile])};

ppparams.perf(1).wcbffile = ['w' ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile];

clear Vcbf