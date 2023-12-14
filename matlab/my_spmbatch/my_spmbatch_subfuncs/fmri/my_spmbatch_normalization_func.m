function [wfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_normalization_func(ne,ppparams,params,delfiles,keepfiles)

Vfunc = spm_vol(ppparams.funcfile{ne});
Rfunc= spm_vol(ppparams.reffunc{ne});

for i=1:numel(Vfunc)
    wfuncfiles{i,1} = [ppparams.funcfile{ne} ',' num2str(i)];
end

%% Normalization of the functional scan
if ne==ppparams.echoes(1)
    funcnormest.subj.vol = {wfuncfiles{1,1}};
    funcnormest.eoptions.biasreg = 0.0001;
    funcnormest.eoptions.biasfwhm = 60;
    funcnormest.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
    funcnormest.eoptions.affreg = 'mni';
    funcnormest.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    funcnormest.eoptions.fwhm = 0;
    funcnormest.eoptions.samp = 3;

    spm_run_norm(funcnormest);

    ppparams.deffile = spm_file(ppparams.funcfile{ne}, 'prefix','y_');
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

defs.comp{1}.def         = {ppparams.deffile};
defs.comp{2}.idbbvox.vox = woptions.vox;
defs.comp{2}.idbbvox.bb  = woptions.bb;
defs.out{1}.pull.fnames  = '';
defs.out{1}.pull.savedir.savesrc = 1;
defs.out{1}.pull.interp  = woptions.interp;
defs.out{1}.pull.mask    = 1;
defs.out{1}.pull.fwhm    = [0 0 0];
defs.out{1}.pull.prefix  = woptions.prefix;

defs.out{1}.pull.fnames = wfuncfiles(:,1);

Nii = nifti(defs.comp{1}.def);
vx  = sqrt(sum(Nii.mat(1:3,1:3).^2));
if det(Nii.mat(1:3,1:3))<0, vx(1) = -vx(1); end

o   = Nii.mat\[0 0 0 1]';
o   = o(1:3)';
dm  = size(Nii.dat);
bb  = [-vx.*(o-1) ; vx.*(dm(1:3)-o)];

defs.comp{2}.idbbvox.vox = woptions.vox;
defs.comp{2}.idbbvox.bb  = woptions.bb;
defs.comp{2}.idbbvox.vox(~isfinite(defs.comp{2}.idbbvox.vox)) = vx(~isfinite(defs.comp{2}.idbbvox.vox));
defs.comp{2}.idbbvox.bb(~isfinite(defs.comp{2}.idbbvox.bb)) = bb(~isfinite(defs.comp{2}.idbbvox.bb));

[~,wfuncdat] = my_spmbatch_deformations(defs);

ppparams.funcfile{ne} = spm_file(ppparams.funcfile{ne}, 'prefix','w');
keepfiles{numel(keepfiles)+1} = {ppparams.funcfile{ne}};