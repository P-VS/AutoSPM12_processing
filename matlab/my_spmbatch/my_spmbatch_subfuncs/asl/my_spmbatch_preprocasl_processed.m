function [delfiles,keepfiles] = my_spmbatch_preprocasl_processed(ppparams,params,delfiles,keepfiles)

if ~params.multiple_scans, params.num_scans = [1]; end

for i=params.num_scans
    if params.multiple_scans
        ppparams.subm0scan = fullfile(ppparams.subasldir ,[ppparams.substring '_asl-' num2str(i) '_m0scan.nii']);
        ppparams.subm0json = fullfile(ppparams.subasldir ,[ppparams.substring '_asl-' num2str(i) '_m0scan.json']);
        ppparams.subdeltam = fullfile(ppparams.subasldir ,[ppparams.substring '_asl-' num2str(i) '_deltam.nii']);
        ppparams.subdmjson = fullfile(ppparams.subasldir ,[ppparams.substring '_asl-' num2str(i) '_deltam.json']);

        if params.aslge.cbfmap_present
            ppparams.subcbf = fullfile(ppparams.subasldir ,[ppparams.substring '_asl-' num2str(i) '_cbf.nii']);
        end
    else
        ppparams.subm0scan = fullfile(ppparams.subasldir ,[ppparams.substring '_asl_m0scan.nii']);
        ppparams.subm0json = fullfile(ppparams.subasldir ,[ppparams.substring '_asl_m0scan.json']);
        ppparams.subdeltam = fullfile(ppparams.subasldir ,[ppparams.substring '_asl_deltam.nii']);
        ppparams.subdmjson = fullfile(ppparams.subasldir ,[ppparams.substring '_asl_deltam.json']);

        if params.aslge.cbfmap_present
            ppparams.subcbf = fullfile(ppparams.subasldir ,[ppparams.substring '_asl_cbf.nii']);
        end
    end

    if params.reorient
        [pth nm ext] = fileparts(ppparams.subm0scan);
        transfile = fullfile(ppparams.subasldir,[nm '_reorient.mat']);
        if isfile(transfile)
            load(transfile,'M')
            transM = M;
        else        
            transM = my_spmbatch_vol_set_com(ppparams.subm0scan);
            transM(1:3,4) = -transM(1:3,4);
        end
    
        Vm0 = spm_vol(ppparams.subm0scan);
        MM = Vm0.private.mat0;
    
        Vm0 = my_reset_orientation(Vm0,transM*MM);
    
        m0dat = spm_read_vols(Vm0);

        m0mask = my_spmbatch_mask(m0dat);
        m0dat(m0mask<0.5) = 0;
    
        Vm0.fname = spm_file(ppparams.subm0scan, 'prefix','r');
        Vm0.descrip = 'reoriented';
        Vm0 = spm_create_vol(Vm0);
        Vm0 = spm_write_vol(Vm0,m0dat);
    
        ppparams.subm0scan = spm_file(ppparams.subm0scan, 'prefix','r');
        delfiles{numel(delfiles)+1} = {ppparams.subm0scan};

        auto_acpc_reorient(ppparams.subm0scan,'PD');

        Vm0 = spm_vol(ppparams.subm0scan);
        MM = Vm0.mat;
    
        Vmd = spm_vol(ppparams.subdeltam);
    
        Vmd = my_reset_orientation(Vmd,MM);
    
        mddat = spm_read_vols(Vmd);
    
        Vmd.fname = spm_file(ppparams.subdeltam, 'prefix','r');
        Vmd.descrip = 'reoriented';
        Vmd = spm_create_vol(Vmd);
        Vmd = spm_write_vol(Vmd,mddat);
    
        ppparams.subdeltam = spm_file(ppparams.subdeltam, 'prefix','r');
        delfiles{numel(delfiles)+1} = {ppparams.subdeltam};
    
        if params.aslge.cbfmap_present
            Vcbf = spm_vol(ppparams.subcbf);
        
            Vcbf = my_reset_orientation(Vcbf,MM);
        
            cbfdat = spm_read_vols(Vcbf);
        
            Vcbf.fname = spm_file(ppparams.subcbf, 'prefix','r');
            Vcbf.descrip = 'reoriented';
            Vcbf = spm_create_vol(Vcbf);
            Vcbf = spm_write_vol(Vcbf,cbfdat);
        
            ppparams.subcbf = spm_file(ppparams.subcbf, 'prefix','r');
            delfiles{numel(delfiles)+1} = {ppparams.subcbf};
        end
    end

    %% Corregistration M0map to deltam
    estwrite.ref(1) = {ppparams.subdeltam};
    estwrite.source(1) = {ppparams.subm0scan};
    estwrite.other = {ppparams.subm0scan};
    estwrite.eoptions.cost_fun = 'nmi';
    estwrite.eoptions.sep = [4 2];
    estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    estwrite.eoptions.fwhm = [7 7];
    estwrite.roptions.interp = 4;
    estwrite.roptions.wrap = [0 0 0];
    estwrite.roptions.mask = 0;
    estwrite.roptions.prefix = 'r';
    
    out_coreg = spm_run_coreg(estwrite);

    ppparams.subm0scan = spm_file(ppparams.subm0scan, 'prefix','r');
    delfiles{numel(delfiles)+1} = {ppparams.subm0scan};

    %% T1 correction of M0
    if params.aslge.do_cbfmapping
        [ppparams,delfiles,keepfiles] = my_spmbatch_asl_M0correction(ppparams,params,delfiles,keepfiles);
    end

    %% Quantification of CBF
    if params.aslge.do_cbfmapping
        [ppparams,delfiles,keepfiles] = my_spmbatch_asl_cbfquantification(ppparams,delfiles,keepfiles);
    end

    %% Normalization of the ASL scan
    if params.asl.do_normalization
        if params.aslge.cbfmap_present || params.aslge.do_cbfmapping
            resamplescans = {ppparams.subm0scan,ppparams.subdeltam,ppparams.subcbf};
        else
            resamplescans = {ppparams.subm0scan,ppparams.subdeltam};
        end

        aslnormestwrite.subj.vol = {ppparams.subm0scan};
        aslnormestwrite.subj.resample = resamplescans;
        aslnormestwrite.eoptions.biasreg = 0.0001;
        aslnormestwrite.eoptions.biasfwhm = 60;
        aslnormestwrite.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
        aslnormestwrite.eoptions.affreg = 'mni';
        aslnormestwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        aslnormestwrite.eoptions.fwhm = 0;
        aslnormestwrite.eoptions.samp = 3;
        aslnormestwrite.woptions.bb = [-78 -112 -70;78 76 85];
        aslnormestwrite.woptions.vox = params.asl.normvox;
        aslnormestwrite.woptions.interp = 1;
        aslnormestwrite.woptions.prefix = 'w';
    
        spm_run_norm(aslnormestwrite);
    
        ppparams.deffile = spm_file(ppparams.subm0scan, 'prefix','y_');
        delfiles{numel(delfiles)+1} = {ppparams.deffile};

        keepfiles{numel(keepfiles)+1} = {spm_file(ppparams.subm0scan, 'prefix','w')};
        keepfiles{numel(keepfiles)+1} = {spm_file(ppparams.subdeltam, 'prefix','w')};
    
        if params.aslge.cbfmap_present || params.aslge.do_cbfmapping
            keepfiles{numel(keepfiles)+1} = {spm_file(ppparams.subcbf, 'prefix','w')};
        end
    end
end
