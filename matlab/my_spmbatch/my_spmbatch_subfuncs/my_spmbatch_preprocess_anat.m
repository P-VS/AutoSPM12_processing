function [subanat,c1im,c2im,c3im,deffile,delfiles,keepfiles] = my_spmbatch_preprocess_anat(substring,subanatdir,preproc_anat,params,delfiles,keepfiles)

if ~exist(fullfile(subanatdir,[substring '_T1w_Crop_1.nii']),'file')
    nsubannat = [substring '_T1w.nii'];
    nsubanstring = [substring '_T1w'];
else
    nsubannat = [substring '_T1w_Crop_1.nii'];
    nsubanstring = [substring '_T1w_Crop_1'];
end
subanat = fullfile(subanatdir,nsubannat);

%% Segmentation and normalization of the anatomical T1w scan

if exist(fullfile(preproc_anat,['wr' nsubannat]),'file')
    if exist(fullfile(preproc_anat,['wc1r' nsubannat]),'file')
        if exist(fullfile(preproc_anat,['wc2r' nsubannat]),'file')
            if exist(fullfile(preproc_anat,['wc3r' nsubannat]),'file')
                params.do_segmentation = 0;

                c1im = fullfile(preproc_anat,['wc1r' nsubannat]);
                c2im = fullfile(preproc_anat,['wc2r' nsubannat]);
                c3im = fullfile(preproc_anat,['wc3r' nsubannat]);
            end
        end
    end
end
if exist(fullfile(preproc_anat,['wr' nsubanstring]),'file')
    if exist(fullfile(preproc_anat,['w' nsubanstring '_brain_pve_0.nii']),'file')
        if exist(fullfile(preproc_anat,['w' nsubanstring '_brain_pve_1.nii']),'file')
            if exist(fullfile(preproc_anat,['w' nsubanstring '_brain_pve_2.nii']),'file')
                params.do_segmentation = 0;

                c1im = fullfile(preproc_anat,['w' nsubanstring '_brain_pve_1.nii']);
                c2im = fullfile(preproc_anat,['w' nsubanstring '_brain_pve_2.nii']);
                c3im = fullfile(preproc_anat,['w' nsubanstring '_brain_pve_0.nii']);
            end
        end
    end
end

if params.reorient
    [pth nm ext] = fileparts(subanat);
    transfile = fullfile(subanatdir,[nm '_reorient.mat']);
    if isfile(transfile)
        load(transfile,'M')
        transM = M;
    else
        transM = eye(4);
    end

    Vanat = spm_vol(subanat);
    MM = Vanat.private.mat0;

    Vanat = my_reset_orientation(Vanat,transM * MM);

    anatdat = spm_read_vols(Vanat);

    Vanat.fname = spm_file(subanat, 'prefix','r');
    Vanat.descrip = 'reoriented';
    Vanat = spm_create_vol(Vanat);
    Vanat = spm_write_vol(Vanat,anatdat);

    subanat = spm_file(subanat, 'prefix','r');
    delfiles{numel(delfiles)+1} = {subanat};

    auto_acpc_reorient(subanat,'T1');
end

if params.do_segmentation
    %%Do segmentation

    preproc.channel.vols = {subanat};
    preproc.channel.biasreg = 0.001;
    preproc.channel.biasfwhm = 60;
    preproc.channel.write = [0 0];
    preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')};
    preproc.tissue(1).ngaus = 1;
    preproc.tissue(1).native = [1 0];
    preproc.tissue(1).warped = [0 0];
    preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
    preproc.tissue(2).ngaus = 1;
    preproc.tissue(2).native = [1 0];
    preproc.tissue(2).warped = [1 0];
    preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
    preproc.tissue(3).ngaus = 2;
    preproc.tissue(3).native = [1 0];
    preproc.tissue(3).warped = [0 0];
    preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
    preproc.tissue(4).ngaus = 3;
    preproc.tissue(4).native = [0 0];
    preproc.tissue(4).warped = [0 0];
    preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
    preproc.tissue(5).ngaus = 4;
    preproc.tissue(5).native = [0 0];
    preproc.tissue(5).warped = [0 0];
    preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
    preproc.tissue(6).ngaus = 2;
    preproc.tissue(6).native = [0 0];
    preproc.tissue(6).warped = [0 0];
    preproc.warp.mrf = 1;
    preproc.warp.cleanup = 1;
    preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    preproc.warp.affreg = 'mni';
    preproc.warp.fwhm = 0;
    preproc.warp.samp = 3;
    preproc.warp.write = [0 1];
    preproc.warp.vox = NaN;
    preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];

    spm_preproc_run(preproc);

    c1im = spm_file(subanat, 'prefix','c1');
    c2im = spm_file(subanat, 'prefix','c2');
    c3im = spm_file(subanat, 'prefix','c3');

    if params.do_normalization
        %%Normalization of the T1w anatomical scan and the segmentation maps
        deffile = spm_file(subanat, 'prefix','y_');

        segnormwrite.subj.def(1) = {deffile};
        segnormwrite.subj.resample = {subanat,c1im,c2im,c3im};
        segnormwrite.woptions.bb = [-78 -112 -70;78 76 85];
        segnormwrite.woptions.vox = params.normvox;
        segnormwrite.woptions.interp = 4;
        segnormwrite.woptions.prefix = 'w';

        spm_run_norm(segnormwrite);

        delfiles{numel(delfiles)+1} = {deffile};
        delfiles{numel(delfiles)+1} = {c1im};
        delfiles{numel(delfiles)+1} = {c2im};
        delfiles{numel(delfiles)+1} = {c3im};

        keepfiles{numel(keepfiles)+1} = {spm_file(subanat, 'prefix','w')};
        keepfiles{numel(keepfiles)+1} = {spm_file(c1im, 'prefix','w')};
        keepfiles{numel(keepfiles)+1} = {spm_file(c2im, 'prefix','w')};
        keepfiles{numel(keepfiles)+1} = {spm_file(c3im, 'prefix','w')};

        c1im = spm_file(subanat, 'prefix','wc1');
        c2im = spm_file(subanat, 'prefix','wc2');
        c3im = spm_file(subanat, 'prefix','wc3');
    else
        keepfiles{numel(keepfiles)+1} = {c1im};
        keepfiles{numel(keepfiles)+1} = {c2im};
        keepfiles{numel(keepfiles)+1} = {c3im};
    end
elseif params.do_normalization
    %%Normalization of the T1w anatomical scan

    anatnormestwrite.subj.vol = {subanat};
    anatnormestwrite.subj.resample = {subanat};
    anatnormestwrite.eoptions.biasreg = 0.0001;
    anatnormestwrite.eoptions.biasfwhm = 60;
    anatnormestwrite.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
    anatnormestwrite.eoptions.affreg = 'mni';
    anatnormestwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    anatnormestwrite.eoptions.fwhm = 0;
    anatnormestwrite.eoptions.samp = 3;
    anatnormestwrite.woptions.bb = [-78 -112 -70;78 76 85];
    anatnormestwrite.woptions.vox = params.normvox;
    anatnormestwrite.woptions.interp = 4;
    anatnormestwrite.woptions.prefix = 'w';

    spm_run_norm(anatnormestwrite);

    deffile = spm_file(subanat, 'prefix','y_');
    delfiles{numel(delfiles)+1} = {deffile};
    keepfiles{numel(keepfiles)+1} = {spm_file(subanat, 'prefix','w')};
end

if ~exist('deffile','var'); deffile=''; end