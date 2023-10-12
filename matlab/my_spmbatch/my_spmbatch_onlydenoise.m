function [delfiles,keepfiles] = my_spmbatch_onlydenoise(sub,ses,task,datpath,params)

delfiles = {};
keepfiles = {};

ppparams = params;

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.task = task;

ppparams.substring = ['sub-' num2str(sub,'%02d')];

if ~isfolder(fullfile(datpath,ppparams.substring))
    ppparams.substring = ['sub-' num2str(sub,'%03d')];
end

ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%03d')]);

if ~isfolder(ppparams.subpath)
    ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%02d')]);
end

ppparams.preproc_anat = fullfile(ppparams.subpath,'preproc_anat');

%% Get functional data

if ppparams.mecombined, ppparams.echoes = [1]; end

ppparams.subfmridir = fullfile(ppparams.subpath,'func');

ppparams.subfuncstring = [ppparams.substring '_task-' ppparams.task '_bold'];

if ~params.meepi
    ppparams.funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring '.json']);
else
    ppparams.funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring '_e1.json']);
end

ppparams.ppfuncdir = fullfile(ppparams.subpath,params.save_folder);

ppparams.rp_file = fullfile(ppparams.ppfuncdir,['rp_e' ppparams.subfuncstring '.txt']);

ppparams.subfuncstring = [params.prefix ppparams.subfuncstring];

if params.do_ICA_AROMA || params.do_DENN
    %% Make masks
    
    if ~ppparams.mecombined
        nsubfuncstring = ['s' ppparams.subfuncstring '_e' num2str(ppparams.echoes(1))];
    else
        nsubfuncstring = ['s' ppparams.subfuncstring];
    end
    
    [Vmask,fmaskdat,~] = my_spmbatch_readSEfMRI(nsubfuncstring,ppparams.ppfuncdir,0,params,10);
    
    % Functional mask
    func_mask = my_spmbatch_mask(fmaskdat);
    
    Vfuncmask = Vmask(1);
    Vfuncmask.fname = fullfile(ppparams.ppfuncdir,['fmask_' nsubfuncstring '.nii']); % Change name to contain subjectID
    Vfuncmask.descrip = 'funcmask';
    Vfuncmask = rmfield(Vfuncmask, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
    spm_write_vol(Vfuncmask, func_mask);
    
    ppparams.fmask = Vfuncmask.fname;
    
    keepfiles{numel(keepfiles)+1} = {Vfuncmask.fname};  
end

%% Do segmentation of func data

if params.do_aCompCor || params.do_ICA_AROMA || params.do_DENN
    if ~ppparams.mecombined
        nsubfuncstring = [ppparams.subfuncstring '_e1'];
    else
        nsubfuncstring = ppparams.subfuncstring;
    end
    
    segfuncfile = fullfile(ppparams.ppfuncdir,[nsubfuncstring '.nii,1']);
    
    preproc.channel.vols = {segfuncfile};
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
    
    segfuncfile = fullfile(ppparams.ppfuncdir,[nsubfuncstring '.nii']);
    
    ppparams.wc1im = spm_file(segfuncfile, 'prefix','c1');
    ppparams.wc2im = spm_file(segfuncfile, 'prefix','c2');
    ppparams.wc3im = spm_file(segfuncfile, 'prefix','c3');
    
    keepfiles{numel(keepfiles)+1} = {ppparams.wc1im};  
    keepfiles{numel(keepfiles)+1} = {ppparams.wc2im}; 
    keepfiles{numel(keepfiles)+1} = {ppparams.wc3im}; 
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir,[nsubfuncstring '._seg8.mat'])};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir,['y_' nsubfuncstring '.nii'])};    
end

%% Denoising with motion derivatives and squared regressors (24 regressors)
if params.do_mot_derivatives
    fprintf('Start motion derivatives \n')

    confounds = load(ppparams.rp_file);
    confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
    confounds = cat(2,confounds,power(confounds,2));

    ppparams.rp_file = spm_file(ppparams.rp_file, 'prefix','der_','ext','.txt');

    writematrix(confounds,ppparams.rp_file,'Delimiter','tab');

    keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
end

%% Denoising with aCompCor covariates
if params.do_aCompCor
    
    for ie=ppparams.echoes
        if ~ppparams.mecombined
            nsubfuncstring = [ppparams.subfuncstring '_e' num2str(ie)];
        else
            nsubfuncstring = ppparams.subfuncstring;
        end
        
        [tedata{ie}.Vfunc,tedata{ie}.wfuncdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(nsubfuncstring,ppparams.ppfuncdir,0,params,Inf);
    end

    fprintf('Start aCompCor \n')

    if ~ppparams.mecombined
        nechoes = numel(ppparams.echoes);

        for i=1:nechoes
            if i==1
                voldim = tedata{ie}.Vfunc.dim;
                tefuncdat = zeros(voldim(1),voldim(2),voldim(3),numel(tedata{ie}.Vfunc),nechoes);
            end
        
            tefuncdat(:,:,:,:,i) = tedata{ppparams.echoes(i)}.wfuncdat;
        end

        wfuncdat = sum(tefuncdat,5) ./ nechoes;
    else
        wfuncdat = tedata{1}.wfuncdat;
    end

    [ppparams,keepfiles] = my_spmbatch_acompcor(wfuncdat,ppparams,params,keepfiles);

    fprintf('Done aCompCor \n')
end  

%% Denoising with ICA-AROMA
if params.do_ICA_AROMA

    for ie=ppparams.echoes
        if ~ppparams.mecombined
            nsubfuncstring = ['s' ppparams.subfuncstring '_e' num2str(ie)];
        else
            nsubfuncstring = ['s' ppparams.subfuncstring];
        end
        
        [tedata{ie}.Vfunc,tedata{ie}.wfuncdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(nsubfuncstring,ppparams.ppfuncdir,0,params,Inf);
    end

    fprintf('Start ICA-AROMA \n')

    [ppparams,keepfiles,delfiles] = fmri_ica_aroma(tedata,ppparams,keepfiles,delfiles);

    fprintf('Done ICA-AROMA \n')
end

%% Noise regression
if params.do_noiseregression || params.do_bpfilter
    fprintf('Do noise regression \n')

    for ie=ppparams.echoes
        if ~params.do_ICA_AROMA
            if ~ppparams.mecombined
                nsubfuncstring = ['s' ppparams.subfuncstring '_e' num2str(ie)];
            else
                nsubfuncstring = ['s' ppparams.subfuncstring];
            end

            [tedata{ie}.Vfunc,tedata{ie}.wfuncdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(nsubfuncstring,ppparams.ppfuncdir,0,params,Inf);
        end

        [tedata{ie}.wfuncdat,ppparams,delfiles] = my_spmbatch_noiseregression(tedata{ie}.wfuncdat,1,ppparams,params,delfiles);

        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','d');
        keepfiles{numel(keepfiles)+1} = {ppparams.funcfile{ie}};  
    end

    fprintf('Done noise regresion \n')
end