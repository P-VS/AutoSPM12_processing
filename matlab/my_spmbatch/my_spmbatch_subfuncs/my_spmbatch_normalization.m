function [wfuncdat,funcfile,keepfiles] = my_spmbatch_normalization(subanat,reffunc,funcfile,deffile,params,keepfiles)

%% Coregistration func to anat

funccorestimate.ref = {subanat};
funccorestimate.source = {reffunc};
funccorestimate.other = {''};
funccorestimate.eoptions.cost_fun = 'nmi';
funccorestimate.eoptions.sep = [4 2];
funccorestimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
funccorestimate.eoptions.fwhm = [7 7];

%spm_run_coreg(funccorestimate);
x  = spm_coreg(char(funccorestimate.ref), char(funccorestimate.source), funccorestimate.eoptions);
M  = spm_matrix(x);

Vfunc = spm_vol(funcfile);

if ~isempty(Vfunc(1).private.extras) && isstruct(Vfunc(1).private.extras) && isfield(Vfunc(1).private.extras,'mat')
    for i=1:size(Vfunc(1).private.dat,4)
        omat = Vfunc(1).private.extras.mat;
        mat(:,:,i) = M\omat(:,:,i);
    end
end
for k=1:numel(Vfunc)
    MM = Vfunc(k).mat;
    Vfunc(k).mat = M\MM;
    Vfunc(k).private.mat = Vfunc(k).mat;
    if ~isempty(Vfunc(k).private.extras) && isstruct(Vfunc(k).private.extras) && isfield(Vfunc(k).private.extras,'mat')
        Vfunc(k).private.extras.mat = mat;
    end

    Vfunc(k) = spm_create_vol(Vfunc(k));
end

%% Normalise func

for i=1:numel(Vfunc)
    wfuncfiles{i,1} = [funcfile ',' num2str(i)];
end

%Write the spatially normalised data

woptions.bb = [-78 -112 -70;78 76 85];
woptions.vox = params.normvox;
woptions.interp = 4;
woptions.prefix = 'w';

defs.comp{1}.def         = {deffile};
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

funcfile = spm_file(funcfile, 'prefix','w');
keepfiles{numel(keepfiles)+1} = {funcfile};