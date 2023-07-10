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
        if k==1
            Vfunc(k).pinfo = [];
        else
            Vfunc(k).pinfo = Vfunc(1).pinfo;
        end
        Vfunc(k).n = [k 1];
        Vfunc(k) = spm_create_vol(Vfunc(k));
        Vfunc(k) = spm_write_vol(Vfunc(k),funcdat(:,:,:,k));
    end

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
    
    %%
    if params.fieldmap
        [ppparams,delfiles,keepfiles] = my_spmbatch_fieldmap(ie,ppparams,delfiles,keepfiles);
    else
        ppparams.vdm_file = '';
    end

    ppparams.prefix = '';
    %%
    if params.do_realignment 
        [funcdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ie,ppparams,params,keepfiles,delfiles);
        
        Vfunc = spm_vol(ppparams.funcfile{ie});
    end

    %%
    if params.pepolar
        [funcdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(funcdat,numdummy,ie,ppparams,params,delfiles,keepfiles);
        ppparams.prefix = ['u' ppparams.prefix];
    end
    
    if params.do_slicetime    
        %% Slice time correction
    
        fprintf('Do slice time correction\n')
        
        SliceTimes = jsondat.SliceTiming;
        nsl= numel(jsondat.SliceTiming);
        
        funcdat=my_spmbatch_st(funcdat,Vfunc,SliceTimes,tr);
    
        if ~params.meepi || contains(params.combination,'none')
            for k=1:numel(Vfunc)
                Vfunc(k).fname = spm_file(ppparams.funcfile{ie}, 'prefix',['a' ppparams.prefix]);
                Vfunc(k).descrip = 'my_spmbatch - slice time correction';
                if k==1
                    Vfunc(k).pinfo = [];
                else
                    Vfunc(k).pinfo = Vfunc(1).pinfo;
                end
                Vfunc(k).n = [k 1];
                Vfunc(k) = spm_create_vol(Vfunc(k));
                Vfunc(k) = spm_write_vol(Vfunc(k),funcdat(:,:,:,k));
            end
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
                if k==1
                    Vfunc(k).pinfo = [];
                else
                    Vfunc(k).pinfo = Vfunc(1).pinfo;
                end
                Vfunc(k).n = [k 1];
                Vfunc(k) = spm_create_vol(Vfunc(k));
                Vfunc(k) = spm_write_vol(Vfunc(k),funcdat(:,:,:,k));
            end
        end
    
        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix',ppparams.prefix);
        delfiles{numel(delfiles)+1} = {ppparams.funcfile{ie}};
    end
    
    if params.meepi
        tefuncdata{ie}.data = funcdat;
        tefuncdata{ie}.Vfunc = Vfunc;
    end
end
    
%%
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
    %%
    if params.do_normalization
        [wfuncdat,ppparams,keepfiles] = my_spmbatch_normalization(ie,ppparams,params,keepfiles);
    end
       
    %%
    %%Denoising with motion derivatives and squared regressors (24 regressors)
    if params.do_mot_derivatives
        fprintf('Start motion derivatives \n')
    
        confounds = load(ppparams.rp_file);
        confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));
        confounds = cat(2,confounds,power(confounds,2));
    
        ppparams.rp_file = spm_file(ppparams.rp_file, 'prefix','der_','ext','.txt');
    
        writematrix(confounds,ppparams.rp_file,'Delimiter','tab');
    
        keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
    end
 
    %%Denoising with aCompCor covariates
    if params.do_aCompCor  && ie==ppparams.echoes(1)
        fprintf('Start aCompCor \n')
    
        [ppparams,keepfiles] = my_spmbatch_acompcor(wfuncdat,ppparams,params,keepfiles);
    end  

    if params.do_noiseregression || params.do_bpfilter
        fprintf('Do noise regression \n')
    
        [wfuncdat,ppparams,delfiles] = my_spmbatch_noiseregression(wfuncdat,ie,ppparams,params,delfiles);
    end
    
    if params.do_smoothing
        fprintf('Start smoothing \n')
    
        %% Smooth func
        
        if ~exist('wfuncdat','var')
    
            for i=1:numel(Vfunc)
                sfuncfiles{i,1} = [ppparams.funcfile{ie} ',' num2str(i)];
            end
    
            smooth.data = sfuncfiles(:,1);
            smooth.fwhm = [params.smoothfwhm params.smoothfwhm params.smoothfwhm];
            smooth.dtype = 0;
            smooth.im = 0;
            smooth.prefix = 's';
            
            spm_run_smooth(smooth)
        else
            Vfunc = spm_vol(ppparams.funcfile{ie});
    
            sfuncdat = zeros([Vfunc(1).dim(1),Vfunc(1).dim(2),Vfunc(1).dim(3),numel(Vfunc)]);
    
            spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');
    
            for i=1:numel(Vfunc)
                [pth,nm,~] = fileparts(Vfunc(i).fname);
    
                Q = fullfile(pth, ['s' nm  '.nii,' num2str(Vfunc(i).n)]);
                my_spmbatch_smooth(wfuncdat(:,:,:,i),Vfunc(i),Q,[params.smoothfwhm params.smoothfwhm params.smoothfwhm],0);
    
                spm_progress_bar('Set',i);
            end
    
            spm_progress_bar('Clear');
        end
    
        ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','s');
        keepfiles{numel(keepfiles)+1} = {ppparams.funcfile{ie}};    
    
        fprintf('Done smoothing \n')
    end
end