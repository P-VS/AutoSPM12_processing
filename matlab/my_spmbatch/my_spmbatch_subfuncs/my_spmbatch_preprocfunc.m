function [delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles)

if ~params.meepi
    ppparams.subfuncstring = [ppparams.substring '_task-' ppparams.task '_bold'];
    params.echoes = [1];
else
    ppparams.subfuncstring = [ppparams.substring '_task-' ppparams.task '_bold_e'];
end

if ~params.meepi
    ppparams.funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring '.json']);
else
    ppparams.funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring '1.json']);
end

jsondat = fileread(ppparams.funcjsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;

numdummy = floor(params.dummytime/tr);

if ppparams.testscript
    readvols = 20;  %Only for a quick test of the batch script
else
    readvols = Inf;
end
ppparams.echoes = params.echoes;
ppparams.meepi = params.meepi;

for ie=ppparams.echoes

    %% Load func data without dummy scans

    if ppparams.meepi
        nsubfuncstring = [ppparams.subfuncstring num2str(ie)];
    else
        nsubfuncstring = ppparams.subfuncstring;
    end

    [Vfunc,funcdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(nsubfuncstring,ppparams.subfmridir,numdummy,params,readvols);

    fprintf('Write data without dummys\n')
    
    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(ppparams.funcfile{ie}, 'prefix','e');
        Vfunc(k).descrip = 'my_spmbatch - remove dummys';
        Vfunc(k).n = [k 1];
    end

    Vfunc = myspm_write_vol_4d(Vfunc,funcdat);

    ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','e');
    delfiles{numel(delfiles)+1} = {ppparams.funcfile{ie}};
    
    if params.reorient   
        if ie==ppparams.echoes(1)
            auto_acpc_reorient([ppparams.funcfile{ie} ',1'],'EPI');
    
            Vfunc = spm_vol(ppparams.funcfile{ie});
            MM = Vfunc(1).mat;
        end
    
        Vfunc = my_reset_orientation(Vfunc,MM);
        for k=1:numel(Vfunc)
            Vfunc(k) = spm_create_vol(Vfunc(k));
        end
    end
    
    Rfunc = Vfunc(1);
    rdat = spm_read_vols(Rfunc);
    
    Rfunc.fname = spm_file(ppparams.funcfile{ie}, 'prefix','f');
    Rfunc.descrip = 'my_spmbatch - first volume';
    Rfunc.n = [1 1];
    Rfunc = spm_create_vol(Rfunc);
    Rfunc = spm_write_vol(Rfunc,rdat);
    
    ppparams.reffunc{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','f');
    delfiles{numel(delfiles)+1} = {ppparams.reffunc{ie}};
    
    %% Fieldmap geometric correction
    if params.fieldmap
        [ppparams,delfiles,keepfiles] = my_spmbatch_fieldmap(ie,ppparams,delfiles,keepfiles);
    else
        ppparams.vdm_file = '';
    end

    ppparams.prefix = '';

    %% Realignment
    if params.do_realignment 
        [funcdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ie,ppparams,params,keepfiles,delfiles);
        
        Vfunc = spm_vol(ppparams.funcfile{ie});
    end

    %% Topup geometric correction
    if params.pepolar
        [funcdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(funcdat,numdummy,ie,ppparams,params,delfiles,keepfiles);
        ppparams.prefix = ['u' ppparams.prefix];
    end
    
    %% Slice time correction
    if params.do_slicetime    
    
        fprintf('Do slice time correction\n')
        
        SliceTimes = jsondat.SliceTiming;
        nsl= numel(jsondat.SliceTiming);
        
        funcdat=my_spmbatch_st(funcdat,Vfunc,SliceTimes,tr);
    
        if ~params.meepi || contains(params.combination,'none')
            for k=1:numel(Vfunc)
                Vfunc(k).fname = spm_file(ppparams.funcfile{ie}, 'prefix',['a' ppparams.prefix]);
                Vfunc(k).descrip = 'my_spmbatch - slice time correction';
                Vfunc(k).n = [k 1];
            end

            Vfunc = myspm_write_vol_4d(Vfunc,funcdat);
        else
            ppparams.prefix = ['a' ppparams.prefix];
        end
    
        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix',['a' ppparams.prefix]);
        delfiles{numel(delfiles)+1} = {ppparams.funcfile{ie}};
    else
        if ~params.meepi || contains(params.combination,'none')
            for k=1:numel(Vfunc)
                Vfunc(k).fname = spm_file(ppparams.funcfile{ie}, 'prefix',ppparams.prefix);
                Vfunc(k).descrip = 'my_spmbatch';
                Vfunc(k).n = [k 1];
            end

            Vfunc = myspm_write_vol_4d(Vfunc,funcdat);
        end
    
        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix',ppparams.prefix);
        delfiles{numel(delfiles)+1} = {ppparams.funcfile{ie}};
    end
    
    if params.meepi
        tefuncdata{ie}.data = funcdat;
        tefuncdata{ie}.Vfunc = Vfunc;
    end
end
    
%% Combine multiple TE timeseries for ME-fMRI
if params.meepi & ~contains(params.combination,'none')
    fprintf('Combine echoes\n')

    [~,Vfunc] = my_spmbatch_combineMEfMRI(tefuncdata,ppparams,params);

    ppparams.funcfile{1} = Vfunc(1).fname;
    delfiles{numel(delfiles)+1} = {ppparams.funcfile{1}};

    Rfunc = Vfunc(1);
    rdat = spm_read_vols(Rfunc);
    
    Rfunc.fname = spm_file(ppparams.funcfile{1}, 'prefix','f');
    Rfunc.descrip = 'my_spmbatch - first volume';
    Rfunc.n = [1 1];
    Rfunc = spm_create_vol(Rfunc);
    Rfunc = spm_write_vol(Rfunc,rdat);
    
    ppparams.reffunc{1} = spm_file(ppparams.funcfile{1}, 'prefix','f');
    delfiles{numel(delfiles)+1} = {ppparams.reffunc{1}};

    ppparams.echoes = [1];
    ppparams.meepi = false;
end  

for ie=ppparams.echoes

    %% Normalization of func data
    if params.do_normalization
        [wfuncdat,ppparams,keepfiles] = my_spmbatch_normalization(ie,ppparams,params,keepfiles);
    end
     
    if ie==ppparams.echoes(1)

        %% Extend motion regressos with derivatives and squared regressors (24 regressors) for denoising
        if params.do_mot_derivatives
            fprintf('Start motion derivatives \n')
        
            confounds = load(ppparams.rp_file);
            confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
            confounds = cat(2,confounds,power(confounds,2));
        
            ppparams.rp_file = spm_file(ppparams.rp_file, 'prefix','der_','ext','.txt');
        
            writematrix(confounds,ppparams.rp_file,'Delimiter','tab');
        
            keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
        end
             
        %% Do segmentation of func data       
        if params.do_aCompCor || params.do_ICA_AROMA || params.do_DENN
               
            [sfpath,sfname,~] = fileparts(ppparams.funcfile{ie});
           
            segfuncfile = fullfile(sfpath,[sfname '.nii,1']);
        
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
                
            ppparams.wc1im = spm_file(ppparams.funcfile{ie}, 'prefix','c1');
            ppparams.wc2im = spm_file(ppparams.funcfile{ie}, 'prefix','c2');
            ppparams.wc3im = spm_file(ppparams.funcfile{ie}, 'prefix','c3');
    
            keepfiles{numel(keepfiles)+1} = {ppparams.wc1im};  
            keepfiles{numel(keepfiles)+1} = {ppparams.wc2im}; 
            keepfiles{numel(keepfiles)+1} = {ppparams.wc3im}; 
            delfiles{numel(delfiles)+1} = {fullfile(sfpath,[sfname '._seg8.mat'])};
            delfiles{numel(delfiles)+1} = {fullfile(sfpath,['y_' sfname '.nii'])};
        end
        
        %% Denoising with aCompCor covariates
        if params.do_aCompCor
            
            fprintf('Start aCompCor \n')
        
            [ppparams,keepfiles] = my_spmbatch_acompcor(wfuncdat,ppparams,params,keepfiles);
        
            fprintf('Done aCompCor \n')
        end  

    end 

    %% Smooth func        
    if params.do_smoothing
        fprintf('Start smoothing \n')
    
        if ~exist('wfuncdat','var')
            [~,sfname,~] = fileparts(ppparams.funcfile{ie});

            [Vfunc,wfuncdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(sfname,ppparams.subfmridir,0,params,Inf);
    
        end

        Vfunc = spm_vol(ppparams.funcfile{ie});

        sfuncdat = zeros([Vfunc(1).dim(1),Vfunc(1).dim(2),Vfunc(1).dim(3),numel(Vfunc)]);

        spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');

        for i=1:numel(Vfunc)
            [pth,nm,~] = fileparts(Vfunc(i).fname);

            Q = fullfile(pth, ['s' nm  '.nii,' num2str(Vfunc(i).n)]);
            swfuncdat = my_spmbatch_smooth(wfuncdat(:,:,:,i),Vfunc(i),Q,[params.smoothfwhm params.smoothfwhm params.smoothfwhm],0);

            sfuncdat(:,:,:,i) = swfuncdat;

            clear swfuncdat

            spm_progress_bar('Set',i);
        end

        spm_progress_bar('Clear');

        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','s');
        keepfiles{numel(keepfiles)+1} = {ppparams.funcfile{ie}};    

        tedata{ie}.Vfunc = spm_vol(ppparams.funcfile{ie});
        tedata{ie}.wfuncdat = sfuncdat;

        fprintf('Done smoothing \n')
    end
    
end

if params.do_ICA_AROMA || params.do_DENN
    %% Make masks
    
    % Functional mask
    fmaskdat = tedata{ppparams.echoes(1)}.wfuncdat;
    
    func_mask = my_spmbatch_mask(fmaskdat);
    
    Vfuncmask = tedata{ppparams.echoes(1)}.Vfunc(1);
    Vfuncmask.fname = spm_file(ppparams.funcfile{ppparams.echoes(1)}, 'prefix','fmask'); 
    Vfuncmask.descrip = 'funcmask';
    Vfuncmask = rmfield(Vfuncmask, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
    spm_write_vol(Vfuncmask, func_mask);
    
    ppparams.fmask = Vfuncmask.fname;
    ppparams.ppfuncdir = ppparams.subfmridir;
    
    keepfiles{numel(keepfiles)+1} = {Vfuncmask.fname};  
end

%% Denoising with ICA-AROMA
if params.do_ICA_AROMA
    fprintf('Start ICA-AROMA \n')

    [ppparams,keepfiles,delfiles] = fmri_ica_aroma(tedata,ppparams,keepfiles,delfiles);

    fprintf('Done ICA-AROMA \n')
end

%% Noise regression
if params.do_noiseregression || params.do_bpfilter
    fprintf('Do noise regression \n')

    for ie=ppparams.echoes
        
        [tedata{ie}.wfuncdat,ppparams,delfiles] = my_spmbatch_noiseregression(tedata{ie}.wfuncdat,1,ppparams,params,delfiles);

        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','d');
        keepfiles{numel(keepfiles)+1} = {ppparams.funcfile{ie}};  
    end

    fprintf('Done noise regresion \n')
end
