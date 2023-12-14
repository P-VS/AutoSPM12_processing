function [delfiles,keepfiles] = my_spmbatch_normfslsegments_job(sub,ses,datpath,params)

substring = ['sub-' num2str(sub,'%02d')];

if ~isfolder(fullfile(datpath,substring))
    substring = ['sub-' num2str(sub,'%03d')];
end

subpath = fullfile(datpath,substring,['ses-' num2str(ses,'%03d')]);

if ~isfolder(subpath)
    subpath = fullfile(datpath,substring,['ses-' num2str(ses,'%02d')]);
end

subanatdir = fullfile(subpath,'anat');

preproc_anat = fullfile(subpath,'preproc_anat');

delfiles = {};
keepfiles = {};

if ~exist(fullfile(subanatdir,[substring '_T1w_Crop_1.nii']),'file')
    nsubannat = [substring '_T1w.nii'];
    nm = [substring '_T1w'];
else
    nsubannat = [substring '_T1w_Crop_1.nii'];
    nm = [substring '_T1w_Crop_1'];
end
subanat = fullfile(subanatdir,nsubannat);

gmim = fullfile(subanatdir,[nm '_brain_pve_1.nii']);
wmim = fullfile(subanatdir,[nm '_brain_pve_2.nii']);
csfim = fullfile(subanatdir,[nm '_brain_pve_0.nii']);

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {subanat};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {gmim;wmim;csfim};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = params.normvox;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

deffile = spm_file(subanat, 'prefix','y_');
delfiles{numel(delfiles)+1} = {deffile};

%Run matlabbatch
spm_jobman('run', matlabbatch);

keepfiles{numel(keepfiles)+1} = {spm_file(gmim, 'prefix','w')};
keepfiles{numel(keepfiles)+1} = {spm_file(wmim, 'prefix','w')};
keepfiles{numel(keepfiles)+1} = {spm_file(csfim, 'prefix','w')};

end

function P = my_reset_orientation(P,MM)

    if ~isempty(P(1).private.extras) && isstruct(P(1).private.extras) && isfield(P(1).private.extras,'mat')
        for i=1:size(P(1).private.dat,4)
            mat(:,:,i) = MM;
        end
    end
    for k=1:numel(P)
        P(k).mat = MM;
        P(k).private.mat = MM;
        if ~isempty(P(k).private.extras) && isstruct(P(k).private.extras) && isfield(P(k).private.extras,'mat')
            P(k).private.extras.mat = mat;
        end
    end

end