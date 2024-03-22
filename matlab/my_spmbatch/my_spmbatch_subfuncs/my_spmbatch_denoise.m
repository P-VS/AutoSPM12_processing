function [delfiles,keepfiles] = my_spmbatch_denoise(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

ppparams = params.denoise;

if contains(ppparams.prefix,'c'), ppparams.mecombined = true; else ppparams.mecombined = false; end

if ~ppparams.meepi
    ppparams.echoes = [1];
    ppparams.mecombined = true;
elseif ppparams.mecombined 
    ppparams.echoes = [1];
end

ppparams.reorient = false;

if contains(ppparams.prefix(1),'s') && numel(ppparams.prefix)>1
    ppparams.prefix = ppparams.prefix(2:end);
end

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.task = task;
ppparams.use_parallel = params.use_parallel;
ppparams.save_intermediate_results = params.save_intermediate_results;

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring)), ppparams.substring = ['sub-' num2str(sub,'%03d')]; end

ppparams.sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring,ppparams.sesstring)), ppparams.sesstring = ['ses-' num2str(ses,'%03d')]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.subfuncdir = fullfile(ppparams.subpath,'func');

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.ppfuncdir = fullfile(ppparams.subpath,params.save_folder);

if ~isfolder(ppparams.ppfuncdir)
    fprintf(['No preprocessed func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

%% Get original functional data

jnamefilters(1).name = ppparams.substring;
jnamefilters(1).required = true;

jnamefilters(2).name = ppparams.sesstring;
jnamefilters(2).required = false;

jnamefilters(3).name = ['run-' num2str(ppparams.run)];
if params.func.mruns, jnamefilters(3).required = true; else jnamefilters(3).required = false; end

jnamefilters(4).name = ['task-' ppparams.task];
jnamefilters(4).required = true;

jnamefilters(5).name = '_bold';
jnamefilters(5).required = true;

if ppparams.meepi
    jnamefilters(6).name = '_echo-1';
    jnamefilters(6).required = true;
end

funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);

if isempty(funcjsonlist) && ppparams.meepi
    jnamefilters(6).name = '_e1';
    jnamefilters(6).required = true;
    
    funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
end

if isempty(funcjsonlist)
    fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

ppparams.funcjsonfile = fullfile(funcjsonlist(1).folder,funcjsonlist(1).name);

%% Get preprocessed functional data

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(ppparams.run)];
if params.denoise.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['task-' ppparams.task];
namefilters(4).required = true;

namefilters(5).name = '_bold';
namefilters(5).required = true;

funcniilist = my_spmbatch_dirfilelist(ppparams.ppfuncdir,'nii',namefilters,false);

if isempty(funcniilist)
    fprintf(['No nii files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

for ie=ppparams.echoes
    if ppparams.meepi && ~ppparams.mecombined %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(ie)]),contains({funcniilist.name},['_e' num2str(ie)])));
        if isempty(tmp)
            fprintf(['no fmri data found for echo ' num2str(ie) '\n'])
            fprintf('\nPP_Error\n');
            return
        end
    
        edirniilist = funcniilist(tmp);
    else
        edirniilist = funcniilist;
    end

    prefixlist = split({edirniilist.name},'sub-');
    if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

    tmp = find(strcmp(prefixlist,ppparams.prefix));
    if ~isempty(tmp), ppparams.func(ie).wfuncfile = edirniilist(tmp).name; end

    tmp = find(strcmp(prefixlist,['s' ppparams.prefix]));
    if ~isempty(tmp), ppparams.func(ie).sfuncfile = edirniilist(tmp).name; end

    tmp = find(strcmp(prefixlist,['fmask_s' ppparams.prefix]));
    if ~isempty(tmp), ppparams.fmask = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

    if ie==ppparams.echoes(1)
        tmp = find(strcmp(prefixlist,['c1' ppparams.prefix]));
        if ~isempty(tmp), ppparams.wc1im = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

        tmp = find(strcmp(prefixlist,['c2' ppparams.prefix]));
        if ~isempty(tmp), ppparams.wc2im = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end

        tmp = find(strcmp(prefixlist,['c3' ppparams.prefix]));
        if ~isempty(tmp), ppparams.wc3im = fullfile(edirniilist(tmp).folder,edirniilist(tmp).name); end
    end

    %% Smooth func 
    if ~isfield(ppparams.func(ie),'sfuncfile')
        if ~isfield(ppparams.func(ie),'wfuncfile')
            fprintf(['no preprocessed fmri data found for echo ' num2str(ie) '\n'])
            fprintf('\nPP_Error\n');
            return
        end

        fprintf('Do smoothing \n')
    
        Vfunc = spm_vol(fullfile(ppparams.ppfuncdir,ppparams.func(ie).wfuncfile));
        wfuncdat = spm_read_vols(Vfunc);

        spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');

        for i=1:numel(Vfunc)
            [pth,nm,~] = fileparts(Vfunc(i).fname);

            Q = fullfile(pth, ['s' nm  '.nii,' num2str(Vfunc(i).n)]);
            my_spmbatch_smooth(wfuncdat(:,:,:,i),Vfunc(i),Q,[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

            spm_progress_bar('Set',i);
        end

        spm_progress_bar('Clear');

        ppparams.func(ie).sfuncfile = ['s' ppparams.func(ie).wfuncfile];   

        clear wfuncdat

        fprintf('Done smoothing \n')
    end
end

%% Get rp_... file
rnamefilters(1).name = ppparams.substring;
rnamefilters(1).required = true;

rnamefilters(2).name = ppparams.sesstring;
rnamefilters(2).required = false;

rnamefilters(3).name = ['run-' num2str(ppparams.run)];
if params.func.mruns, rnamefilters(3).required = true; else rnamefilters(3).required = false; end

rnamefilters(4).name = ['task-' ppparams.task];
rnamefilters(4).required = true;

rnamefilters(5).name = '_bold';
rnamefilters(5).required = true;

rnamefilters(6).name = 'rp_';
rnamefilters(6).required = true;

funcrplist = my_spmbatch_dirfilelist(ppparams.ppfuncdir,'txt',rnamefilters,false);

if isempty(funcrplist)
    fprintf(['No rp_... files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

ppparams.rp_file = fullfile(funcrplist(1).folder,funcrplist(1).name);

%% Get der_... file
if params.denoise.do_mot_derivatives
    rnamefilters(6).name = 'der_';
    rnamefilters(6).required = true;
    
    funcderlist = my_spmbatch_dirfilelist(ppparams.ppfuncdir,'txt',rnamefilters,false);
    
    if ~isempty(funcderlist)
        ppparams.der_file = fullfile(funcderlist(1).folder,funcderlist(1).name);
    end
end

%% Get acc_... file
if params.denoise.do_aCompCor 
    rnamefilters(6).name = 'acc_';
    rnamefilters(6).required = true;
    
    funcacclist = my_spmbatch_dirfilelist(ppparams.ppfuncdir,'txt',rnamefilters,false);
    
    if ~isempty(funcacclist)
        ppparams.acc_file = fullfile(funcacclist(1).folder,funcacclist(1).name);
    end
end

%% Get ica_... file
if params.denoise.do_ICA_AROMA 
    rnamefilters(6).name = 'ica_';
    rnamefilters(6).required = true;
    
    funcicalist = my_spmbatch_dirfilelist(ppparams.ppfuncdir,'txt',rnamefilters,false);
    
    if ~isempty(funcicalist)
        ppparams.ica_file = fullfile(funcicalist(1).folder,funcicalist(1).name);
    end
end

%% Make masks
if params.denoise.do_ICA_AROMA || params.denoise.do_DENN
    if ~isfield(ppparams,'fmask')
        fprintf('Make mask \n')
    
        [Vmask,fmaskdat] = my_spmbatch_readSEfMRI(ppparams.ppfuncdir,ppparams.func(1).sfuncfile,0,ppparams,10);
        
        % Functional mask
        func_mask = my_spmbatch_mask(fmaskdat);
        
        Vfuncmask = Vmask(1);
        Vfuncmask.fname = fullfile(ppparams.ppfuncdir,['fmask_' ppparams.func(1).sfuncfile]); % Change name to contain subjectID
        Vfuncmask.descrip = 'funcmask';
        Vfuncmask = rmfield(Vfuncmask, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
        spm_write_vol(Vfuncmask, func_mask);
        
        ppparams.fmask = Vfuncmask.fname;
        
        keepfiles{numel(keepfiles)+1} = {Vfuncmask.fname};  
    end
end

%% Do segmentation of func data
if params.denoise.do_aCompCor || params.denoise.do_ICA_AROMA || params.denoise.do_DENN
    if ~isfield(ppparams.func(1),'wfuncfile')
        fprintf(['no preprocessed fmri data found for echo ' num2str(1) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    if ~isfield(ppparams,'wc1im') || ~isfield(ppparams,'wc2im') || ~isfield(ppparams,'wc3im')
        fprintf('Do segmentation \n')
        
        segfuncfile = fullfile(ppparams.ppfuncdir,[ppparams.func(1).wfuncfile ',1']);
        
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
        
        segfuncfile = fullfile(ppparams.ppfuncdir,ppparams.func(1).wfuncfile);
        
        ppparams.wc1im = spm_file(segfuncfile, 'prefix','c1');
        ppparams.wc2im = spm_file(segfuncfile, 'prefix','c2');
        ppparams.wc3im = spm_file(segfuncfile, 'prefix','c3');
        
        keepfiles{numel(keepfiles)+1} = {ppparams.wc1im};  
        keepfiles{numel(keepfiles)+1} = {ppparams.wc2im}; 
        keepfiles{numel(keepfiles)+1} = {ppparams.wc3im}; 
    
        sname = split(ppparams.func(1).wfuncfile,'.nii');
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir,[sname{1} '._seg8.mat'])};
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.ppfuncdir,['y_' sname{1} '.nii'])};  
    end
end

%% Denoising with motion derivatives and squared regressors (24 regressors)
if params.denoise.do_mot_derivatives && ~isfield(ppparams,'der_file')
    fprintf('Start motion derivatives \n')

    confounds = load(ppparams.rp_file);
    confounds = confounds(:,1:6);
    confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
    confounds = cat(2,confounds,power(confounds,2));

    ppparams.der_file = spm_file(ppparams.rp_file, 'prefix','der_','ext','.txt');

    writematrix(confounds,ppparams.der_file,'Delimiter','tab');

    keepfiles{numel(keepfiles)+1} = {ppparams.der_file};
end

%% Denoising with aCompCor covariates
if params.denoise.do_aCompCor && ~isfield(ppparams,'acc_file')
    for ie=ppparams.echoes
        [tedata{ie}.Vfunc,tedata{ie}.wfuncdat] = my_spmbatch_readSEfMRI(ppparams.ppfuncdir,ppparams.func(ie).wfuncfile,0,ppparams,Inf);
    end

    fprintf('Start aCompCor \n')
    
    if ppparams.meepi && ~ppparams.mecombined
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
if params.denoise.do_ICA_AROMA && ~isfield(ppparams,'ica_file')
    fprintf('Start ICA-AROMA \n')

    [ppparams,keepfiles,delfiles] = fmri_ica_aroma(ppparams,keepfiles,delfiles);

    fprintf('Done ICA-AROMA \n')
end

%% Noise regression
if params.denoise.do_noiseregression || params.denoise.do_bpfilter
    fprintf('Do noise regression \n')

    for ie=ppparams.echoes
        [~,wfuncdat] = my_spmbatch_readSEfMRI(ppparams.ppfuncdir,ppparams.func(ie).sfuncfile,0,ppparams,Inf);

        [~,ppparams,keepfiles] = my_spmbatch_noiseregression(wfuncdat,ie,ppparams,params,keepfiles);

        clear wfuncdat
    end

    fprintf('Done noise regresion \n')
end