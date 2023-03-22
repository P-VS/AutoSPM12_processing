function [delfiles,keepfiles] = my_spmbatch(sub,ses,task,datpath,params,save_intermediate_results)

testscript = false;

substring = ['sub-' num2str(sub,'%02d')];

if ~isfolder(fullfile(datpath,substring))
    substring = ['sub-' num2str(sub,'%03d')];
end

subpath = fullfile(datpath,substring,['ses-' num2str(ses,'%03d')]);

if ~isfolder(subpath)
    subpath = fullfile(datpath,substring,['ses-' num2str(ses,'%02d')]);
end

subanatdir = fullfile(subpath,'anat');
subfmridir = fullfile(subpath,'func');

preproc_anat = fullfile(subpath,'preproc_anat');

delfiles = {};
keepfiles = {};

if ~exist(fullfile(subanatdir,[substring '_T1w_Crop_1.nii']),'file')
    nsubannat = [substring '_T1w.nii'];
    nsubanstring = [substring '_T1w'];
else
    nsubannat = [substring '_T1w_Crop_1.nii'];
    nsubanstring = [substring '_T1w_Crop_1'];
end
subanat = fullfile(subanatdir,nsubannat);

if params.do_aCompCor
    params.do_segmentation = 1;
end

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

if params.nechoes==1 && ~params.do_onlydenoise
    %% Loading the fMRI time series and deleting dummy scans
    fprintf('Reading the data \n')

    funcjsonfile = fullfile(subfmridir,[substring '_task-' task '_bold.json']);
    funcfile = fullfile(subfmridir,[substring '_task-' task '_bold.nii']);
    
    jsondat = fileread(funcjsonfile);
    jsondat = jsondecode(jsondat);
    tr = jsondat.RepetitionTime;
    
    numdummy = floor(params.dummytime/tr);
    
    Vfunc = spm_vol(funcfile);

    if params.reorient
        transfile = fullfile(subfmridir,[substring '_task-' task '_bold_reorient.mat']);
        if isfile(transfile)
            load(transfile,'M')
            transM = M;
        else
            transM = eye(4);
        end

        MM = Vfunc(1).private.mat0;

        Vfunc = my_reset_orientation(Vfunc,transM * MM);
    end

    Vfunc = Vfunc(numdummy+1:end);
    if testscript
        Vfunc = Vfunc(1:50); %Only for a quick test of the batch script
    end

    funcdat = spm_read_vols(Vfunc);
elseif params.nechoes>1 && ~params.do_onlydenoise
    %% Loading the fMRI time series and deleting dummy scans
    fprintf('Reading the data \n')

    for i=1:params.nechoes
        funcjsonfile = fullfile(subfmridir,[substring '_task-' task '_bold_e' num2str(i) '.json']);
        funcfile = fullfile(subfmridir,[substring '_task-' task '_bold_e' num2str(i) '.nii']);

        jsondat = fileread(funcjsonfile);
        jsondat = jsondecode(jsondat);

        te(i) = jsondat.EchoTime;
        tr = jsondat.RepetitionTime;
    
        numdummy = floor(params.dummytime/tr);

        Vfunc = spm_vol(funcfile);

        if params.reorient
            transfile = fullfile(subfmridir,[substring '_task-' task '_bold_e1_reorient.mat']);
            if isfile(transfile)
                load(transfile,'M')
                transM = M;
            else
                transM = eye(4);
            end

            MM = Vfunc(1).private.mat0;

            Vfunc = my_reset_orientation(Vfunc,transM * MM);
        end
        
        Vfunc = Vfunc(numdummy+1:end);
        if testscript
            Vfunc = Vfunc(1:50); %Only for a quick test of the batch script
        end

        if i==1
            voldim = Vfunc.dim;
            funcdat = zeros(voldim(1),voldim(2),voldim(3),numel(Vfunc));
        end
    
        efuncdat = spm_read_vols(Vfunc);
        funcdat = funcdat+efuncdat.*te(i);
    end

    funcdat = funcdat./sum(te);

    fname=Vfunc(1).fname;
    nfname = split(fname,['_e' num2str(params.nechoes)]);

    for j=1:numel(Vfunc)
        Vfunc(j).fname = [nfname{1} '.nii'];
    end

    funcfile = [nfname{1} '.nii'];
end

if params.do_slicetime    
    %% Slice time correction

    fprintf('Do slice time correction\n')
    
    SliceTimes = jsondat.SliceTiming;
    nsl= numel(jsondat.SliceTiming);
    
    funcdat=my_spmbatch_st(funcdat,Vfunc,SliceTimes,tr);

    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(funcfile, 'prefix','ae');
        Vfunc(k).descrip = 'my_spmbatch - slice time correction';
        Vfunc(k).n = [k 1];
        Vfunc(k) = spm_create_vol(Vfunc(k));
        Vfunc(k) = spm_write_vol(Vfunc(k),funcdat(:,:,:,k));
    end

    funcfile = spm_file(funcfile, 'prefix','ae');
    delfiles{numel(delfiles)+1} = {funcfile};
elseif ~params.do_onlydenoise
    fprintf('Write data without dummys\n')
    
    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(funcfile, 'prefix','e');
        Vfunc(k).descrip = 'my_spmbatch - remove dummys';
        Vfunc(k).n = [k 1];
        Vfunc(k) = spm_create_vol(Vfunc(k));
        Vfunc(k) = spm_write_vol(Vfunc(k),funcdat(:,:,:,k));
    end

    funcfile = spm_file(funcfile, 'prefix','e');
    delfiles{numel(delfiles)+1} = {funcfile};
end

if params.reorient    
    
    auto_acpc_reorient([funcfile ',1'],'EPI');

    Vfunc = spm_vol(funcfile);
    MM = Vfunc(1).mat;

    Vfunc = my_reset_orientation(Vfunc,MM);
    for k=1:numel(Vfunc)
        Vfunc(k) = spm_create_vol(Vfunc(k));
    end
end

if ~params.do_onlydenoise
    Rfunc = Vfunc(1);
    rdat = spm_read_vols(Rfunc);
    
    Rfunc.fname = spm_file(funcfile, 'prefix','f');
    Rfunc.descrip = 'my_spmbatch - first volume';
    Rfunc.n = [1 1];
    Rfunc = spm_create_vol(Rfunc);
    Rfunc = spm_write_vol(Rfunc,rdat);
    
    reffunc = spm_file(funcfile, 'prefix','f');
    delfiles{numel(delfiles)+1} = {reffunc};
end


mbstep = 1;

%%
if params.pepolar
    if params.nechoes==1
        ppfunc = fullfile(subpath,'fmap',[substring '_dir-pi_epi.nii']);

        Vppfunc = spm_vol(ppfunc);
    
        if params.reorient
            transfile = fullfile(subpath,'fmap',[substring '_dir-pi_epi_reorient.mat']);
            if isfile(transfile)
                load(transfile,'M')
                transM = M;
            else
                transM = eye(4);
            end
    
            MM = Vppfunc(1).private.mat0;
            Vppfunc = my_reset_orientation(Vppfunc,transM * MM);
        end
        
        Vppfunc = Vppfunc(numdummy+1);

        ppfuncdat = spm_read_vols(Vppfunc);
    else
        for i=1:params.nechoes
            ppfunc = fullfile(subpath,'fmap',[substring '_dir-pi_epi_e' num2str(i) '.nii']);
    
            Vppfunc = spm_vol(ppfunc);

            if params.reorient
                transfile = fullfile(subpath,'fmap',[substring '_dir-pi_epi_e1_reorient.mat']);
                if isfile(transfile)
                    load(transfile,'M')
                    transM = M;
                else
                    transM = eye(4);
                end
        
                MM = Vppfunc(1).private.mat0;
                Vppfunc = my_reset_orientation(Vppfunc,transM * MM);
            end
            
            Vppfunc = Vppfunc(numdummy+1);
    
            if i==1
                voldim = Vppfunc.dim;
                ppfuncdat = zeros(voldim(1),voldim(2),voldim(3));
            end

            eppfuncdat = spm_read_vols(Vppfunc);
            ppfuncdat = ppfuncdat+eppfuncdat.*te(i);
        end

        ppfuncdat = ppfuncdat./sum(te);

        fname=Vppfunc(1).fname;
        nfname = split(fname,['_e' num2str(params.nechoes)]);
    
        Vppfunc(1).fname = [nfname{1} '.nii'];
    
        ppfunc = [nfname{1} '.nii'];
    end

    Vppfunc.fname = spm_file(ppfunc, 'prefix','f');
    Vppfunc.descrip = 'my_spmbatch - first volume';
    Vppfunc.n = [1 1];
    Vppfunc = spm_create_vol(Vppfunc);
    Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);

    ppfunc = spm_file(ppfunc, 'prefix','f');
    delfiles{numel(delfiles)+1} = {ppfunc};

    if params.reorient
        auto_acpc_reorient(ppfunc,'EPI');
    end

    %% coregister fmap to func

    rfuncpsstep = mbstep;
    
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {ppfunc};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other = {ppfunc};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    mbstep = mbstep+1;

    delfiles{numel(delfiles)+1} = {ppfunc};
    delfiles{numel(delfiles)+1} = {spm_file(ppfunc, 'prefix','r')};

    %% HYSCO fieldmap
    
    pedir = jsondat.PhaseEncodingDirection;
    
    if contains(pedir,'i')
        pedim = 1;
        WrapD = [1 0 0];
    elseif contains(pedir,'j')
        pedim = 2;
        WrapD = [0 1 0];
    else
        pedim = 3;
        WrapD = [0 0 1];
    end
    
    hyscostep = mbstep;
    
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_up(1) = {reffunc};
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_dw(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{rfuncpsstep}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_up(1) = {''};
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_dw(1) = {''};
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.perm_dim = pedim;
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.dummy_fast = 1;
    matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.outdir = {''};
    
    mbstep = mbstep+1;

    delfiles{numel(delfiles)+1} = {spm_file(Vppfunc(1).fname, 'prefix','u2r')};
    delfiles{numel(delfiles)+1} = {spm_file(Vppfunc(1).fname, 'prefix','HySCov2_r')};
    delfiles{numel(delfiles)+1} = {spm_file(reffunc, 'prefix','u2')};

    %% Convert fieldmap into vdm
    
    te = jsondat.EchoTime*1000;
    trt = jsondat.TotalReadoutTime*1000;
    
    if contains(pedir,'-')
        blipdir=-1;
    else
        blipdir=1;
    end
    
    fmstep = mbstep;
    
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = cfg_dep('HySCO: Inhomogeneity field', substruct('.','val', '{}',{hyscostep}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fieldmap'));
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap = {''};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [te te];
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdir;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = trt;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {reffunc};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    
    mbstep = mbstep+1;

    d = dir(fullfile(subpath,'func'));
    ddir=[d(:).isdir];
    dfolders = {d(ddir).name};
    delfolders = find(contains(dfolders,'derivatives'));

    if ~isempty(delfolders)
        for delf=1:numel(delfolders)
            rmdir(fullfile(subpath,'func',dfolders{delfolders(delf)}),'s');
        end
    end

    [pth,name,ext] = fileparts(reffunc);
    vdm_file = fullfile(subpath,'func','derivatives','HySCO-Run',['vdm5_' name '_desc-H-HySCO-ESTIMATED-FIELDMAP_dwi.nii']);

    delfiles{numel(delfiles)+1} = {fullfile(subpath,'func','derivatives')}; 
    delfiles{numel(delfiles)+1} = {spm_file(reffunc, 'prefix','u')};

elseif params.fieldmap
    %%Load fieldmap data per echo and coregister to func
    
    e1dat = fullfile(subpath,'fmap',[substring '_fmap_echo-1_am.nii']);
    e1json = fullfile(subpath,'fmap',[substring '_fmap_echo-1_am.json']);
    
    if isfile(e1dat)
        e1phdat = fullfile(subpath,'fmap',[substring '_fmap_echo-1_ph.nii']);

        e1jsondat = fileread(e1json);
        e1jsondat = jsondecode(e1jsondat);
        te1 = e1jsondat.EchoTime*1000;
    
        Ve1amp = spm_vol(e1dat);
        Ve1ph  = spm_vol(e1phdat);
    else
        e1dat = fullfile(subpath,'fmap',[substring '_fmap_echo-1.nii']);
        e1json = fullfile(subpath,'fmap',[substring '_fmap_echo-1.json']);

        e1jsondat = fileread(e1json);
        e1jsondat = jsondecode(e1jsondat);
        te1 = e1jsondat.EchoTime*1000;
    
        Ve1 = spm_vol(e1dat);
        Ve1amp = Ve1(1);
        Ve1ph  = Ve1(2);
    
        Ve1amp = do_save_intermediate_results(Ve1amp,spm_read_vols(Ve1amp),'am','spm - amplitude');
        Ve1ph = do_save_intermediate_results(Ve1ph,spm_read_vols(Ve1ph),'ph','spm - phase');
    
        delfiles{numel(delfiles)+1} = {Ve1amp(1).fname};
        delfiles{numel(delfiles)+1} = {Ve1ph(1).fname};
    end

    fme1step = mbstep;

    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {Ve1amp(1).fname};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other(1) = {Ve1ph(1).fname};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    mbstep = mbstep+1;

    amp1file = spm_file(Ve1amp(1).fname, 'prefix','r');
    ph1file = spm_file(Ve1ph(1).fname, 'prefix','r');

    delfiles{numel(delfiles)+1} = {amp1file};
    delfiles{numel(delfiles)+1} = {ph1file};

    e2dat = fullfile(subpath,'fmap',[substring '_fmap_echo-2_am.nii']);
    e2json = fullfile(subpath,'fmap',[substring '_fmap_echo-2_am.json']);

    if isfile(e2dat)
        e2phdat = fullfile(subpath,'fmap',[substring '_fmap_echo-2_ph.nii']);

        e2jsondat = fileread(e2json);
        e2jsondat = jsondecode(e2jsondat);
        te2 = e2jsondat.EchoTime*1000;
    
        Ve2amp = spm_vol(e2dat);
        Ve2ph  = spm_vol(e2phdat);
    else
        e2dat = fullfile(subpath,'fmap',[substring '_fmap_echo-2.nii']);
        e2json = fullfile(subpath,'fmap',[substring '_fmap_echo-2.json']);

        e2jsondat = fileread(e2json);
        e2jsondat = jsondecode(e2jsondat);
        te2 = e2jsondat.EchoTime*1000;
    
        Ve2 = spm_vol(e2dat);
        Ve2amp = Ve2(1);
        Ve2ph  = Ve2(2);
    
        Ve2amp = do_save_intermediate_results(Ve2amp,spm_read_vols(Ve2amp),'am','spm - amplitude');
        Ve2ph = do_save_intermediate_results(Ve2ph,spm_read_vols(Ve2ph),'ph','spm - phase');
    
        delfiles{numel(delfiles)+1} = {Ve2amp(1).fname};
        delfiles{numel(delfiles)+1} = {Ve2ph(1).fname};
    end

    fme2step = mbstep;

    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {Ve2amp(1).fname};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other(1) = {Ve2ph(1).fname};
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    mbstep = mbstep+1;

    amp2file = spm_file(Ve2amp(1).fname, 'prefix','r');
    ph2file = spm_file(Ve2ph(1).fname, 'prefix','r');

    delfiles{numel(delfiles)+1} = {amp2file};
    delfiles{numel(delfiles)+1} = {ph2file};

    %% Fieldmap
    
    fmstep = mbstep;
    
    pedir = jsondat.PhaseEncodingDirection;
    
    if contains(pedir,'-')
        blipdim = -1;
    else
        blipdim = 1;
    end
    
    try
        trt = jsondat.TotalReadoutTime*1000;
    catch
        trt = jsondat.AcquisitionMatrixPE*jsondat.EffectiveEchoSpacing*1000;
    end
    
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase(1) = {amp1file};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag(1) = {ph1file};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase(1) = {amp2file};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag(1) = {ph2file};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [te1 te2];
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdim;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = trt;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.session.epi(1) = {reffunc};
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    
    mbstep = mbstep+1;

    vdm_file = spm_file(amp1file, 'prefix','vdm5_sc');

    delfiles{numel(delfiles)+1} = {spm_file(ph1file, 'prefix','m')};
    delfiles{numel(delfiles)+1} = {spm_file(ph1file, 'prefix','bmask')};
    delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','sc')};
    delfiles{numel(delfiles)+1} = {spm_file(amp2file, 'prefix','sc')};
    delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','fpm_sc')};
    delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','vdm5_sc')};
    delfiles{numel(delfiles)+1} = {spm_file(reffunc, 'prefix','u')};
end

%%Run matlabbatch
if exist("matlabbatch",'var')
    spm_jobman('run', matlabbatch);
end

%%
if params.do_realignment

    if params.pepolar || params.fieldmap
        %% Realign and unwarp the func series
        
        realignunwarp.data.scans = {funcfile};
        realignunwarp.data.pmscan(1) = {vdm_file};
        realignunwarp.eoptions.quality = 0.9;
        realignunwarp.eoptions.sep = 4;
        realignunwarp.eoptions.fwhm = 5;
        realignunwarp.eoptions.rtm = 0;
        realignunwarp.eoptions.einterp = 2;
        realignunwarp.eoptions.ewrap = [0 0 0];
        realignunwarp.eoptions.weight = '';
        realignunwarp.uweoptions.basfcn = [12 12];
        realignunwarp.uweoptions.regorder = 1;
        realignunwarp.uweoptions.lambda = 100000;
        realignunwarp.uweoptions.jm = 0;
        realignunwarp.uweoptions.fot = [4 5];
        realignunwarp.uweoptions.sot = [];
        realignunwarp.uweoptions.uwfwhm = 4;
        realignunwarp.uweoptions.rem = 1;
        realignunwarp.uweoptions.noi = 5;
        realignunwarp.uweoptions.expround = 'Average';
        realignunwarp.uwroptions.uwwhich = [2 1];
        realignunwarp.uwroptions.rinterp = 4;
        realignunwarp.uwroptions.wrap = WrapD;
        realignunwarp.uwroptions.mask = 1;
        realignunwarp.uwroptions.prefix = 'u';
        
        spm_run_realignunwarp(realignunwarp);

        rp_file = spm_file(funcfile, 'prefix','rp_','ext','.txt');
    
        keepfiles{numel(keepfiles)+1} = {spm_file(funcfile, 'prefix','rp_','ext','.txt')};
    
        reffunc = spm_file(funcfile, 'prefix','meanu');
        funcfile = spm_file(funcfile, 'prefix','u');
    
        delfiles{numel(delfiles)+1} = {reffunc};
        delfiles{numel(delfiles)+1} = {funcfile};
    else
        %% Reslice the func series
        
        realignestwrite.data{1} = {funcfile};
        realignestwrite.eoptions.quality = 0.9;
        realignestwrite.eoptions.sep = 4;
        realignestwrite.eoptions.fwhm = 5;
        realignestwrite.eoptions.rtm = 1;
        realignestwrite.eoptions.interp = 2;
        realignestwrite.eoptions.wrap = [0 0 0];
        realignestwrite.eoptions.weight = '';
        realignestwrite.roptions.which = [2 1];
        realignestwrite.roptions.interp = 4;
        realignestwrite.roptions.wrap = [0 0 0];
        realignestwrite.roptions.mask = 1;
        realignestwrite.roptions.prefix = 'r';
        
        spm_run_realign(realignestwrite);

        rp_file = spm_file(funcfile, 'prefix','rp_','ext','.txt');
    
        keepfiles{numel(keepfiles)+1} = {spm_file(funcfile, 'prefix','rp_','ext','.txt')};
    
        reffunc = spm_file(funcfile, 'prefix','mean');
        funcfile = spm_file(funcfile, 'prefix','r');
    
        delfiles{numel(delfiles)+1} = {reffunc};
        delfiles{numel(delfiles)+1} = {funcfile};
    end
end

%%
if params.do_normalization

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

end

if params.do_onlydenoise
    fprintf('Start loading preprocessed data \n')

    funcfile = fullfile(subpath,'preproc_func',['wuae' substring '_task-' task '_bold.nii']);

    if ~exist(funcfile)
        funcfile = fullfile(subpath,'preproc_func',['wrae' substring '_task-' task '_bold.nii']);
    end

    Vfunc = spm_vol(funcfile);

    wfuncdat = spm_read_vols(Vfunc);

    rp_file = fullfile(subpath,'preproc_func',['rp_ae' substring '_task-' task '_bold.txt']);

    fprintf('Done loading preprocessed data \n')
end

%%Denoising with motion derivatives and squared regressors (24 regressors)
if params.do_mot_derivatives & exist('rp_file','var')
    fprintf('Start motion derivatives \n')

    if exist(rp_file,'file')
        confounds = load(rp_file);
        confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
        confounds = cat(2,confounds,power(confounds,2));

        rp_file = spm_file(rp_file, 'prefix','der_','ext','.txt');

        writematrix(confounds,rp_file,'Delimiter','tab');

        keepfiles{numel(keepfiles)+1} = {rp_file};
    end

    fprintf('Done motion derivatives \n')
end

%%Denoising with aCompCor covariates
if params.do_aCompCor & exist('rp_file','var')
    fprintf('Start aCompCor \n')

    if exist(rp_file,'file')
        confounds = load(rp_file);
    else
        confounds = [];
    end

    GM = spm_vol(c1im);
    WM = spm_vol(c2im);
    CSF = spm_vol(c3im);

    gmdat = spm_read_vols(GM);
    wmdat = spm_read_vols(WM);
    csfdat = spm_read_vols(CSF);

    braindat = gmdat+wmdat;
    braindat(braindat<0.2)=0;
    braindat(braindat>0.0)=1;

    csfdat(braindat>0.0)=0;
    csfdat(csfdat<0.8)=0;
    csfdat(csfdat>0.0)=1;

    if ~exist('wfuncdat','var')
        acc_confounds = fmri_acompcor(funcfile,{csfdat},params.Ncomponents,'confounds',confounds,'PolOrder',1);
    else
        acc_confounds = fmri_acompcor(wfuncdat(:,:,:,:),{csfdat},params.Ncomponents,'confounds',confounds,'PolOrder',1);
    end

    if exist(rp_file)
        confounds = cat(2,confounds,acc_confounds);
    
        rp_file = spm_file(rp_file, 'prefix','acc_','ext','.txt');
    
        writematrix(confounds,rp_file,'Delimiter','tab');
    else
        rp_file = spm_file(funcfile, 'prefix','acc_','ext','.txt');
    
        writematrix(acc_confounds,rp_file,'Delimiter','tab');
    end

    keepfiles{numel(keepfiles)+1} = {rp_file};

    fprintf('Done aCompCor \n')
end

if params.do_noiseregression || params.do_bpfilter
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
end

if params.do_smoothing
    fprintf('Start smoothing \n')

    %% Smooth func
    
    if ~exist('wfuncdat','var')

        for i=1:numel(Vfunc)
            sfuncfiles{i,1} = [funcfile ',' num2str(i)];
        end

        swarfuncstep = mbstep;
        
        smooth.data = sfuncfiles(:,1);
        smooth.fwhm = [params.smoothfwhm params.smoothfwhm params.smoothfwhm];
        smooth.dtype = 0;
        smooth.im = 0;
        smooth.prefix = 's';
        
        spm_run_smooth(smooth)
    else
        Vfunc = spm_vol(funcfile);

        sfuncdat = zeros([Vfunc(1).dim(1),Vfunc(1).dim(2),Vfunc(1).dim(3),numel(Vfunc)]);

        spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');

        for i=1:numel(Vfunc)
            [pth,nm,xt] = fileparts(Vfunc(i).fname);

            Q = fullfile(pth, ['s' nm  '.nii,' num2str(Vfunc(i).n)]);
            my_spmbatch_smooth(wfuncdat(:,:,:,i),Vfunc(i),Q,[params.smoothfwhm params.smoothfwhm params.smoothfwhm],0);

            spm_progress_bar('Set',i);
        end

        spm_progress_bar('Clear');
    end

    funcfile = spm_file(funcfile, 'prefix','s');
    keepfiles{numel(keepfiles)+1} = {funcfile};    

    fprintf('Done smoothing \n')
end

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