function [delfiles,keepfiles] = my_spmbatch_preprocfasl(ppparams,params,delfiles,keepfiles)

for ie=params.func.echoes
    %% Load asl data and remove dummy scans
    if ~contains(ppparams.asl(ie).aslprefix,'e')
        jsondat = fileread(ppparams.asl(ie).asljson);
        jsondat = jsondecode(jsondat);
        tr = jsondat.RepetitionTime;
        
        numdummy = floor(params.func.dummytime/tr);
    
        [Vasl,fasldat] = my_spmbatch_readSEfMRI(ppparams.subperfdir,ppparams.asl(ie).aslfile,numdummy,ppparams,Inf);
    
        fprintf('Write data without dummys\n')
    
        if contains(params.asl.isM0scan,'last')
            M0 = Vasl(end);
            Vasl = Vasl(1:end-1);
    
            m0dat = fasldat(:,:,:,end);
            fasldat = fasldat(:,:,:,1:end-1);
      
            fparts = split(ppparams.asl(ie).aslfile,'_asl');
            
            M0.fname = fullfile(ppparams.subperfdir,['e' fparts{1} '_m0scan.nii']);
            M0.descrip = 'my_spmbatch - m0scan';
            M0.n = [1 1];
            M0 = spm_create_vol(M0);
            M0 = spm_write_vol(M0,m0dat);

            ppparams.asl(ie).m0scanfile = [fparts{1} '_m0scan.nii'];
            ppparams.asl(ie).em0scanfile = ['e' fparts{1} '_m0scan.nii'];
            ppparams.asl(ie).m0scanprefix = 'e';

            delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(ie).em0scanfile)};
        end
            
        ppparams.asl(ie).easlfile = ['e' ppparams.asl(ie).aslfile];
        ppparams.asl(ie).aslprefix = ['e'];
    
        for k=1:numel(Vasl)
            Vasl(k).fname = fullfile(ppparams.subperfdir,ppparams.asl(ie).easlfile);
            Vasl(k).descrip = 'my_spmbatch - remove dummys';
            Vasl(k).n = [k 1];
        end
    
        Vasl = myspm_write_vol_4d(Vasl,fasldat);
    
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(ie).easlfile)};
        
        if ppparams.reorient   
            if ie==params.func.echoes(1)
                auto_acpc_reorient([fullfile(ppparams.subperfdir,ppparams.asl(ie).easlfile) ',1'],'EPI');
        
                Vasl = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(ie).easlfile));
                ppparams.MM = Vasl(1).mat;
            end
            
            if ~isfield(ppparams,'MM')
                Vasl = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(params.func.echoes(1)).easlfile));
                ppparams.MM = Vasl(1).mat;
            end
        
            Vasl = my_reset_orientation(Vasl,ppparams.MM);
            for k=1:numel(Vasl)
                Vasl(k) = spm_create_vol(Vasl(k));
            end

            if exist('M0','var')
                M0 = my_reset_orientation(M0,ppparams.MM);

                M0 = spm_create_vol(M0);

                clear M0
            end
        end
    end
            
    ppparams.asl(ie).feaslfile = ['f' ppparams.asl(ie).easlfile];
    
    if ~exist(fullfile(ppparams.subperfdir,ppparams.asl(ie).feaslfile),'file')
        if ~exist('Vasl','var')
            Vasl = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(ie).easlfile));
        end
    
        Rasl = Vasl(1);
        rdat = spm_read_vols(Rasl);
        
        Rasl.fname = fullfile(ppparams.subperfdir,ppparams.asl(ie).feaslfile);
        Rasl.descrip = 'my_spmbatch - first volume';
        Rasl.n = [1 1];
        Rasl = spm_create_vol(Rasl);
        Rasl = spm_write_vol(Rasl,rdat);
    
        ppparams.asl(ie).refasl = ppparams.asl(ie).feaslfile;
        
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(ie).feaslfile)};
    end
      
    %% Realignment
    if (params.func.do_realignment && ~contains(ppparams.asl(ie).aslprefix,'r'))
        fprintf('Do realignment\n')

        ppparams.subfuncdir = ppparams.subperfdir;
 
        ppparams.func(ie).efuncfile = ppparams.asl(ie).easlfile;
        ppparams.func(ie).prefix = ppparams.asl(ie).aslprefix;

        if params.func.meepi && ie>params.func.echoes(1), ppparams.func(params.func.echoes(1)).efuncfile = ppparams.asl(params.func.echoes(1)).easlfile; end
    
        [fasldat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ie,ppparams,params,keepfiles,delfiles);
        
        ppparams.asl(ie).aslprefix = ppparams.func(ie).prefix;
        ppparams.asl(ie).raslfile = ppparams.func(ie).rfuncfile;
        ppparams.asl(ie).meanaslfile = ppparams.func(ie).meanfuncfile;
        ppparams.asl(ie).refasl = ppparams.func(ie).meanfuncfile;

        if ie==1, ppparams.asl(1).rp_file = ppparams.func(1).rp_file; end
    
        Vasl = spm_vol(fullfile(ppparams.subperfdir,[ppparams.asl(ie).aslprefix ppparams.asl(ie).aslfile]));
    end

    %% Coregistration of m0scan to asl
    if ~contains(ppparams.asl(ie).m0scanprefix,'r')
        estwrite.ref(1) = {fullfile(ppparams.subperfdir,ppparams.asl(ie).refasl)};
        estwrite.source(1) = {fullfile(ppparams.subperfdir,ppparams.asl(ie).em0scanfile)};
        estwrite.other = {fullfile(ppparams.subperfdir,ppparams.asl(ie).em0scanfile)};
        estwrite.eoptions.cost_fun = 'nmi';
        estwrite.eoptions.sep = [4 2];
        estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        estwrite.eoptions.fwhm = [7 7];
        estwrite.roptions.interp = 4;
        estwrite.roptions.wrap = [0 0 0];
        estwrite.roptions.mask = 0;
        estwrite.roptions.prefix = 'r';
        
        out_coreg = spm_run_coreg(estwrite);

        ppparams.asl(ie).rm0scanfile = ['r' ppparams.asl(ie).em0scanfile];
        ppparams.asl(ie).m0scanprefix = ['r' ppparams.asl(ie).m0scanprefix];

        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(ie).rm0scanfile)};
    end
    
    if ~exist('fasldat','var') && ...
            ((params.func.pepolar && ~contains(ppparams.asl(ie).aslprefix,'u')) ...
            || (params.asl.do_denoise && ~contains(ppparams.asl(ie).aslprefix,'d')) ...
            || (params.asl.do_removebold && ~contains(ppparams.asl(ie).aslprefix,'c'))...
            || (params.func.do_slicetime && ~contains(ppparams.asl(ie).aslprefix,'a')))
    
        fprintf('Reading the asl data\n')
    
        Vasl=spm_vol(fullfile(ppparams.subperfdir,[ppparams.asl(ie).aslprefix ppparams.asl(ie).aslfile]));
    
        fasldat = spm_read_vols(Vasl);
    end

    Vm0scan = spm_vol(fullfile(ppparams.subperfdir,[ppparams.asl(ie).m0scanprefix ppparams.asl(ie).m0scanfile]));
    m0dat = spm_read_vols(Vm0scan);

    %% Topup geometric correction
    if params.func.pepolar && ~contains(ppparams.asl(ie).aslprefix,'u')
        fprintf('Do pepolar\n')

        tfasldat = cat(4,fasldat,m0dat);
    
        jsondat = fileread(ppparams.func(ie).fmapjsonfile);
        jsondat = jsondecode(jsondat);
        tr = jsondat.RepetitionTime;
        
        numdummy = floor(params.func.dummytime/tr);

        ppparams.func(ie).reffunc = ppparams.asl(ie).refasl;
        ppparams.subfuncdir = ppparams.subperfdir;
        ppparams.func(ie).jsonfile = ppparams.asl(ie).asljson;
        ppparams.func(ie).prefix = ppparams.asl(ie).aslprefix;
        ppparams.func(ie).funcfile = ppparams.asl(ie).aslfile;
    
        [tfasldat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(tfasldat,numdummy,ie,ppparams,params,delfiles,keepfiles);

        fasldat = tfasldat(:,:,:,1:numel(Vasl));
        m0dat = tfasldat(:,:,:,end);

        prefix = ['u' ppparams.asl(ie).aslprefix];
        m0prefix = ['u' ppparams.asl(ie).m0scanprefix];
    else
        prefix = ppparams.asl(ie).aslprefix;
        m0prefix = ppparams.asl(ie).m0scanprefix;
    end

    %% Drift and motion correction
    if  params.asl.do_denoise && ~contains(ppparams.asl(ie).aslprefix,'d')
        fprintf('Do trend and motion filtering')

        if isfield(ppparams.asl(1),'rp_file')
            confounds = load(ppparams.asl(1).rp_file);
            confounds = confounds(:,1:6);
            confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
            confounds = cat(2,confounds,power(confounds,2));
        else
            confounds = [];
        end

        s = size(fasldat);
        fasldat = reshape(fasldat(:,:,:,:),[prod(s(1:end-1)),s(end)]);
    
        [fasldat,~] = fmri_cleaning(fasldat(:,:),params.asl.polort,[],confounds,[],'restoremean','on');

        fasldat = reshape(fasldat(:,:),s);

        prefix = ['d' prefix];
    end

    if params.asl.do_removebold && ~contains(ppparams.asl(1).aslprefix,'c')
        tefasldata{ie}.data = fasldat;
        tefasldata{ie}.Vasl = Vasl;
        tefasldata{ie}.prefix = prefix;
    elseif exist('fasldat','var')
        for k=1:numel(Vasl)
            Vasl(k).fname = fullfile(ppparams.subperfdir,[prefix ppparams.asl(ie).aslfile]);
            Vasl(k).descrip = 'my_spmbatch';
            Vasl(k).n = [k 1];
        end

        Vasl = myspm_write_vol_4d(Vasl,fasldat);

        ppparams.asl(ie).aslprefix = prefix;

        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[prefix ppparams.asl(ie).aslfile])};
    end

    if params.asl.do_removebold && ~contains(ppparams.asl(1).m0scanprefix,'c')
        tem0data{ie}.data = m0dat;
        tem0data{ie}.Vasl = Vm0scan;
        tem0data{ie}.prefix = m0prefix;
    elseif exist('m0dat','var')
        Vm0scan.fname = fullfile(ppparams.subperfdir,[m0prefix ppparams.asl(ie).m0scanfile]);
        Vm0scan.descrip = 'my_spmbatch';
        Vm0scan.n = [1 1];

        Vm0scan = myspm_write_vol_4d(Vm0scan,m0dat);

        ppparams.asl(ie).m0scanprefix = m0prefix;

        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[m0prefix ppparams.asl(ie).m0scanfile])};
    end

    if exist('fasldat','var'), clear fasldat; end
    if exist('m0dat','var'), clear m0dat; end
end

%% Romove BOLD from the ASL-BOLD series
if params.asl.do_removebold && ~contains(ppparams.asl(1).aslprefix,'c')
    fprintf('Remove BOLD from ASL data\n')

    [fasldata,Vasl,ppparams,delfiles] = my_spmbatch_aslremovebold(tefasldata,ppparams,params,delfiles,true);

    Rasl = Vasl(1);
    rdat = spm_read_vols(Rasl);
    
    Rasl.fname = spm_file(Vasl(1).fname, 'prefix','f');
    Rasl.descrip = 'my_spmbatch - first volume';
    Rasl.n = [1 1];
    Rasl = spm_create_vol(Rasl);
    Rasl = spm_write_vol(Rasl,rdat);

    ppparams.asl(1).refasl = ['fc' prefix ppparams.asl(1).aslfile];
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(1).refasl)};

    ppparams.echoes = [1];

    ppparams.asl(1).aslprefix = ['c' prefix];

    clear tefasldata
elseif exist('tefasldata','var')
    fasldata = tefasldata{params.func.echoes(1)}.data;
    Vasl = tefasldata{params.func.echoes(1)}.Vasl;

    ppparams.asl(1).caslfile = [ppparams.asl(1).aslprefix ppparams.asl(1).aslfile];
elseif ~isfield(ppparams.asl(1),'caslfile')
    ppparams.asl(1).caslfile = [ppparams.asl(1).aslprefix ppparams.asl(1).aslfile];
end

if params.asl.do_removebold && ~contains(ppparams.asl(1).m0scanprefix,'c')
    fprintf('Remove BOLD from m0scan\n')

    [m0data,Vm0scan,ppparams,delfiles] = my_spmbatch_aslremovebold(tem0data,ppparams,params,delfiles,false);

    ppparams.asl(1).m0scanprefix = ['c' m0prefix];

    clear tem0data
end

%% Make GM, WM masks
if ~isfield(ppparams.asl(1),'c2m0scanfile') || ~isfield(ppparams.asl(1),'c1m0scanfile') || ~isfield(ppparams.asl(1),'c3m0scanfile')
    preproc.channel.vols = {fullfile(ppparams.subperfdir,[ppparams.asl(1).m0scanprefix ppparams.asl(1).m0scanfile])};
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
    preproc.tissue(2).warped = [0 0];
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
    
    ppparams.asl(1).c1m0scanfile = ['c1' ppparams.asl(1).m0scanprefix ppparams.asl(1).m0scanfile];
    ppparams.asl(1).c2m0scanfile = ['c2' ppparams.asl(1).m0scanprefix ppparams.asl(1).m0scanfile];
    ppparams.asl(1).c3m0scanfile = ['c3' ppparams.asl(1).m0scanprefix ppparams.asl(1).m0scanfile];

    fname = split([ppparams.asl(1).m0scanprefix ppparams.asl(1).m0scanfile],'.nii');
     
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[fname{1} '._seg8.mat'])};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,['y_' fname{1} '.nii'])};  
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(1).c1m0scanfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(1).c2m0scanfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(1).c3m0scanfile)};
end

%% Control-Label subtraction to make deltam series
if ~isfield(ppparams,'deltam')
    fprintf('Start subtraction \n')

    if ~exist('fasldata','var')
        fprintf('Reading the asl data\n')
    
        Vasl=spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(1).caslfile));
    
        fasldata = spm_read_vols(Vasl);
    end

    [fasldata,Vasl,ppparams,delfiles] = my_spmbacth_faslsubtraction(fasldata,Vasl,ppparams,params,delfiles);
end

%% Make cbf series
if params.asl.do_cbfmapping || ~isfield(ppparams,'cbfmap')
    if ~isfield(ppparams,'tm0scan')
        fprintf('Start M0 preprocessing \n')
    
      %  [ppparams,delfiles,keepfiles] = my_spmbatch_asl_M0correction(ppparams,params,delfiles,keepfiles);
    end
end