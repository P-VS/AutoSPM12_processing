function [delfiles,keepfiles] = my_spmbatch(sub,ses,task,datpath,params)

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

if params.do_aCompCor
    params.do_segmentation = 1;
end

[subanat,c1im,c2im,c3im,deffile,delfiles,keepfiles] = my_spmbatch_preprocess_anat(substring,subanatdir,preproc_anat,params,delfiles,keepfiles);

if params.nechoes==1 && ~params.do_onlydenoise
    subfuncstring = [substring '_task-' task '_bold'];

    funcjsonfile = fullfile(subfmridir,[subfuncstring '.json']);
    
    jsondat = fileread(funcjsonfile);
    jsondat = jsondecode(jsondat);
    tr = jsondat.RepetitionTime;
    
    numdummy = floor(params.dummytime/tr);

    if testscript
        readvols = 50;
    else
        readvols = Inf;
    end

    [Vfunc,funcdat] = my_spmbatch_readSEfMRI(subfuncstring,subfmridir,numdummy,params,readvols);
elseif params.nechoes>1 && ~params.do_onlydenoise
    subfuncstring = [substring '_task-' task '_bold_e'];

    funcjsonfile = fullfile(subfmridir,[subfuncstring '1.json']);
    
    jsondat = fileread(funcjsonfile);
    jsondat = jsondecode(jsondat);
    tr = jsondat.RepetitionTime;
    
    numdummy = floor(params.dummytime/tr);

    if testscript
        readvols = 50;
    else
        readvols = Inf;
    end

    [Vfunc,funcdat,funcfile] = my_spmbatch_readMEfMRI(subfuncstring,subfmridir,numdummy,params,readvols);
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

%%
if params.pepolar
    [vdm_file,delfiles,keepfiles] = my_spmbatch_pepolar(subpath,substring,task,numdummy,reffunc,params,delfiles,keepfiles);
elseif params.fieldmap
    [vdm_file,delfiles,keepfiles] = my_spmbatch_fieldmap(subpath,substring,task,reffunc,delfiles,keepfiles);
end

%%
if params.do_realignment
    if ~exist('vdm_file','var'); vdm_file=''; end

    [reffunc,funcfile,rp_file,keepfiles,delfiles] = my_spmbatch_realignunwarp(subpath,substring,task,funcfile,vdm_file,params,keepfiles,delfiles);
end

%%
if params.do_normalization
    [wfuncdat,funcfile,keepfiles] = my_spmbatch_normalization(subanat,reffunc,funcfile,deffile,params,keepfiles);
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

%%
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
    if ~exist('wfuncdat','var')
        Vfunc = spm_vol(funcfile);

        wfuncdat = spm_read_vols(Vfunc);
    end

    [rp_file,keepfiles] = my_spmbatch_acompcor(subfmridir,substring,task,wfuncdat,funcfile,rp_file,c1im,c2im,c3im,params,keepfiles);
end

if params.do_noiseregression || params.do_bpfilter
    [wfuncdat,funcfile,keepfiles] = my_spmbatch_noiseregression(subfmridir,substring,task,funcfile,rp_file,wfuncdat,keepfiles);
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