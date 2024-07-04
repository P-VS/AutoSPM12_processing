function [delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles)

for ie=params.func.echoes
    %% Load func data and remove dummy scans
    if ~contains(ppparams.func(ie).prefix,'e')
        jsondat = fileread(ppparams.func(ie).jsonfile);
        jsondat = jsondecode(jsondat);
        tr = jsondat.RepetitionTime;
        
        numdummy = floor(params.func.dummytime/tr);
    
        [Vfunc,funcdat] = my_spmbatch_readSEfMRI(ppparams.subfuncdir,ppparams.func(ie).funcfile,numdummy,ppparams,Inf);
    
        fprintf('Write data without dummys\n')
    
        if params.func.isaslbold
            if contains(params.asl.isM0scan,'last')
                M0 = Vfunc(end);
                Vfunc = Vfunc(1:end-1);
    
                m0dat = funcdat(:,:,:,end);
                funcdat = funcdat(:,:,:,1:end-1);
    
                if ~exist(fullfile(ppparams.subpath,'perf'),'dir'), mkdir(fullfile(ppparams.subpath,'perf')); end
    
                fparts = split(ppparams.func(ie).funcfile,'_bold');
                
                M0.fname = fullfile(ppparams.subpath,'perf',['e' fparts{1} '_m0scan.nii']);
                M0.descrip = 'my_spmbatch - m0scan';
                M0.n = [1 1];
                M0 = spm_create_vol(M0);
                M0 = spm_write_vol(M0,m0dat);
            end
        end
            
        ppparams.func(ie).efuncfile = ['e' ppparams.func(ie).funcfile];
        ppparams.func(ie).prefix = ['e'];

        for k=1:numel(Vfunc)
            Vfunc(k).fname = fullfile(ppparams.subfuncdir,ppparams.func(ie).efuncfile);
            Vfunc(k).descrip = 'my_spmbatch - remove dummys';
            Vfunc(k).n = [k 1];
        end
    
        Vfunc = myspm_write_vol_4d(Vfunc,funcdat);

        clear funcdat

        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ie).efuncfile)};
        
        if ppparams.reorient   
            if ie==params.func.echoes(1)
                auto_acpc_reorient([fullfile(ppparams.subfuncdir,ppparams.func(ie).efuncfile) ',1'],'EPI');
        
                Vfunc = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(ie).efuncfile));
                ppparams.MM = Vfunc(1).mat;
            end
            
            if ~isfield(ppparams,'MM')
                Vfunc = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(params.func.echoes(1)).efuncfile));
                ppparams.MM = Vfunc(1).mat;
            end
        
            Vfunc = my_reset_orientation(Vfunc,ppparams.MM);
            for k=1:numel(Vfunc)
                Vfunc(k) = spm_create_vol(Vfunc(k));
            end
        end
    end
            
    ppparams.func(ie).fefuncfile = ['f' ppparams.func(ie).efuncfile];

    if ~exist(fullfile(ppparams.subfuncdir,ppparams.func(ie).fefuncfile),'file')
        if ~exist('Vfunc','var')
            Vfunc = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(ie).efuncfile));
        end

        Rfunc = Vfunc(1);
        rdat = spm_read_vols(Rfunc);
        
        Rfunc.fname = fullfile(ppparams.subfuncdir,ppparams.func(ie).fefuncfile);
        Rfunc.descrip = 'my_spmbatch - first volume';
        Rfunc.n = [1 1];
        Rfunc = spm_create_vol(Rfunc);
        Rfunc = spm_write_vol(Rfunc,rdat);

        ppparams.func(ie).reffunc = ppparams.func(ie).fefuncfile;

        clear rdat Rfunc Vfunc
        
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ie).fefuncfile)};
    end
        
    %% Fieldmap geometric correction
    if params.func.fieldmap && ~contains(ppparams.func(ie).prefix,'u')
        fprintf('Do field map\n')

        [ppparams,delfiles,keepfiles] = my_spmbatch_fieldmap(ie,ppparams,delfiles,keepfiles);
    else
        ppparams.func(ie).vdm_file = '';
    end

    %% Realignment
    if (params.func.do_realignment && ~contains(ppparams.func(ie).prefix,'r')) || (params.func.fieldmap && ~contains(ppparams.func(ie).prefix,'u'))
        fprintf('Do realignment\n')

        [funcdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ie,ppparams,params,keepfiles,delfiles);
        
        Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
    end

    if ~exist('funcdat','var') && ...
            ((params.func.pepolar && ~contains(ppparams.func(ie).prefix,'u')) ...
            || (params.func.do_slicetime && ~contains(ppparams.func(ie).prefix,'a'))...
            || (params.func.do_echocombination && ~contains(ppparams.func(ie).prefix,'c')))

        fprintf('Reading the functional data\n')

        Vfunc=spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));

        funcdat = spm_read_vols(Vfunc);
    end

    %% Topup geometric correction
    if params.func.pepolar && ~contains(ppparams.func(ie).prefix,'u')
        fprintf('Do pepolar\n')

        jsondat = fileread(ppparams.func(ie).fmapjsonfile);
        jsondat = jsondecode(jsondat);
        tr = jsondat.RepetitionTime;
        
        numdummy = floor(params.func.dummytime/tr);

        [funcdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(funcdat,numdummy,ie,ppparams,params,delfiles,keepfiles);

        prefix = ['u' ppparams.func(ie).prefix];
    else
        prefix = ppparams.func(ie).prefix;
    end
    
    %% Masking of the functional data
    if exist('funcdat','var') && ie==params.func.echoes(1)
        fprintf('Do masking\n')

        funcmask = my_spmbatch_mask(funcdat);
        funcmask(funcmask>0) = 1;
        funcdat = repmat(funcmask,[1,1,1,numel(funcmask(1,1,1,:))]) .* funcdat;

        ppparams.func(ie).fefuncfile = ['f' prefix ppparams.func(ie).funcfile];

        Rfunc = Vfunc(1);
        
        Rfunc.fname = fullfile(ppparams.subfuncdir,ppparams.func(ie).fefuncfile);
        Rfunc.descrip = 'my_spmbatch - first volume';
        Rfunc.n = [1 1];
        Rfunc = spm_create_vol(Rfunc);
        Rfunc = spm_write_vol(Rfunc,funcdat(:,:,:,1));

        ppparams.func(ie).reffunc = Rfunc.fname;
        
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ie).fefuncfile)};
    end

    %% Save realignd series for ASL processing
    if params.func.isaslbold
        fparts = split([ppparams.func(ie).funcfile],'_bold');

        deltamfile = fullfile(ppparams.subpath,'perf',[prefix fparts{1} '_asl.nii']);

        if ~exist(deltamfile,'file')
            if exist('funcdat','var')
                Vasl = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
        
                for it=1:numel(Vasl)
                    Vasl(it).fname = deltamfile;
                    Vasl(it).descrip = 'my_spmbatch - asl';
                    Vasl(it).n = [it 1];
                end
        
                Vasl = myspm_write_vol_4d(Vasl,funcdat);
            else
                copyfile(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]),deltamfile);
            end
        end
    end

    %% Removing the tagging if ASL-BOLD
    if params.func.isaslbold && exist('funcdat','var')
        dim = size(funcdat);

        taglabels = ones(dim(4),1);
        if contains(params.asl.tagorder(2),'con'), taglabels(2:2:dim(4)) = -1; else taglabels(1:2:dim(4)) = -1; end

        s = size(funcdat);
        funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

        [funcdat,~] = fmri_cleaning(funcdat(:,:),0,[],taglabels,[],'restoremean','on');

        funcdat = reshape(funcdat(:,:),s);
    end
    
    %% Slice time correction
    if params.func.do_slicetime  && ~contains(ppparams.func(ie).prefix,'a')
        fprintf('Do slice time correction\n')
            
        jsondat = fileread(ppparams.func(ie).jsonfile);
        jsondat = jsondecode(jsondat);

        tr = jsondat.RepetitionTime;
        nsl=Vfunc(1).dim(3);
        
        if params.func.isaslbold
            if isfield(jsondat,'LabelingDuration'), tr=tr-jsondat.LabelingDuration; else tr=tr-params.asl.LabelingDuration; end
            if isfield(jsondat,'PostLabelDelay'), tr=tr-jsondat.PostLabelDelay; else tr=tr-params.asl.PostLabelDelay; end
        end

        if ~isfield(ppparams,'SliceTimes')
            jsondat = fileread(ppparams.func(ie).jsonfile);
            jsondat = jsondecode(jsondat);
        
            ppparams.tr = jsondat.RepetitionTime;
            nsl=tVfunc(1).dim(3);
            
            if ~params.func.isaslbold && isfield(jsondat,'SliceTiming')
                ppparams.SliceTimes = jsondat.SliceTiming;
            else
                if isfield(jsondat,'MultibandAccelerationFactor')
                    hbf = jsondat.MultibandAccelerationFactor; 
                    nslex = ceil(nsl/hbf);
                    isl = zeros([1,nslex]);
                    isl(1:2:nslex)=[0:1:(nslex-1)/2];
                    isl(2:2:nslex)=[ceil(nslex/2):1:nslex-1];
                    isl=repmat(isl,[1,hbf]);
                else 
                    isl = [1:2:nsl 2:2:nsl];
                    nslex = nsl;
                end
        
                if params.func.isaslbold
                    if isfield(jsondat,'LabelingDuration'), params.asl.LabelingDuration=jsondat.LabelingDuration; end
                    if isfield(jsondat,'PostLabelDelay'), params.asl.PostLabelDelay=jsondat.PostLabelDelay; end
        
                    TA = ppparams.tr-params.asl.LabelingDuration-params.asl.PostLabelDelay;
        
                    ppparams.SliceTimes = params.asl.LabelingDuration+params.asl.PostLabelDelay+isl*TA/nslex;
                else
                    TA = ppparams.tr/nslex;
            
                    ppparams.SliceTimes = isl*TA/nslex;
                end
            end
        end
        
        funcdat=my_spmbatch_st(funcdat,tVfunc,ppparams.SliceTimes,ppparams.tr);

        prefix = ['a' prefix];
    end

    if params.func.do_echocombination && ~contains(ppparams.func(ie).prefix,'c')
        tefuncdata{ie}.data = funcdat;
        tefuncdata{ie}.Vfunc = Vfunc;
    elseif exist('funcdat','var')
        for k=1:numel(Vfunc)
            Vfunc(k).fname = fullfile(ppparams.subfuncdir,[prefix ppparams.func(ie).funcfile]);
            Vfunc(k).descrip = 'my_spmbatch';
            Vfunc(k).n = [k 1];
        end

        Vfunc = myspm_write_vol_4d(Vfunc,funcdat);

        ppparams.func(ie).prefix = prefix;

        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[prefix ppparams.func(ie).funcfile])};
    end

    if exist('funcdat','var'), clear funcdat; end
end
 
%% Combine multiple TE timeseries for ME-fMRI
if params.func.do_echocombination && ~contains(ppparams.func(1).prefix,'c')
    fprintf('Combine echoes\n')

    [~,Vfunc,ppparams,delfiles] = my_spmbatch_combineMEfMRI(tefuncdata,ppparams,params,delfiles);

    Rfunc = Vfunc(1);
    rdat = spm_read_vols(Rfunc);
    
    Rfunc.fname = spm_file(Vfunc(1).fname, 'prefix','f');
    Rfunc.descrip = 'my_spmbatch - first volume';
    Rfunc.n = [1 1];
    Rfunc = spm_create_vol(Rfunc);
    Rfunc = spm_write_vol(Rfunc,rdat);

    ppparams.func(1).reffunc = ['fc' ppparams.func(1).prefix ppparams.func(1).funcfile];
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(1).reffunc)};

    ppparams.echoes = [1];
    ppparams.meepi = false;

    ppparams.func(1).prefix = ['c' ppparams.func(1).prefix];

    clear tefuncdata
elseif params.func.do_echocombination
    Vfunc = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(1).cfuncfile));
    Rfunc = Vfunc(1);
    rdat = spm_read_vols(Rfunc);
    
    Rfunc.fname = spm_file(Vfunc(1).fname, 'prefix','f');
    Rfunc.descrip = 'my_spmbatch - first volume';
    Rfunc.n = [1 1];
    Rfunc = spm_create_vol(Rfunc);
    Rfunc = spm_write_vol(Rfunc,rdat);

    clear rdat

    ppparams.func(1).reffunc = ['f' ppparams.func(1).cfuncfile];
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(1).reffunc)};

    ppparams.echoes = [1];
    ppparams.meepi = false;

    nfname = split(ppparams.func(1).funcfile,'_echo-');
    ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];

    ppparams.func(1).prefix = ['c' ppparams.func(1).prefix];

else
    ppparams.echoes = params.func.echoes;
    ppparams.meepi = params.func.meepi;
end  

for ie=ppparams.echoes

    %% Normalization of func data
    if params.func.do_normalization && ~contains(ppparams.func(ie).prefix,'w')
        fprintf('Do normalization\n')

        funcfile = fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]);

        if isfield(ppparams.func(ie),'reffuncfile')
            reffunc = fullfile(ppparams.subfuncdir,ppparams.func(ie).reffuncfile);
        else
            reffunc = [funcfile ',1'];
        end

        [wfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_oldnormalization_func(ie,funcfile,reffunc,ppparams,params,delfiles,keepfiles);

        ppparams.func(ie).prefix = ['w' ppparams.func(ie).prefix];
    end
     
    %% Smooth func        
    if params.func.do_smoothing && ~contains(ppparams.func(ie).prefix,'s')
        fprintf('Do smoothing \n')
    
        Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));

        Vout = Vfunc;
        swfuncdat = zeros([Vfunc(1).dim(1),Vfunc(1).dim(2),Vfunc(1).dim(3),numel(Vfunc)]);

        spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');

        for i=1:numel(Vfunc)

            if ~exist('wfuncdat','var'), twfuncdat = spm_read_vols(Vfunc(i)); else twfuncdat = wfuncdat(:,:,:,i); end

            tswfuncdat = my_spmbatch_smooth(twfuncdat,Vfunc(i),[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

            swfuncdat(:,:,:,i) = tswfuncdat;

            spm_progress_bar('Set',i);

            clear tswfuncdat twfuncdat
        end

        spm_progress_bar('Clear');
    
        for j=1:numel(Vout)
            Vout(j).fname = fullfile(ppparams.subfuncdir,['s' ppparams.func(ie).prefix ppparams.func(ie).funcfile]);
            Vout(j).descrip = 'my_spmbatch - smooth';
            Vout(j).pinfo = [1,0,0];
            Vout(j).n = [j 1];
        end
    
        Vout = myspm_write_vol_4d(Vout,swfuncdat);

        ppparams.func(ie).sfuncfile = ['s' ppparams.func(ie).prefix ppparams.func(ie).funcfile];
        keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ie).sfuncfile)};    

        ppparams.func(ie).prefix = ['s' ppparams.func(ie).prefix];

        fprintf('Done smoothing \n')
    end
    
    if ~exist('wfuncdat','var'), clear wfuncdat; end
end

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile])}; 
