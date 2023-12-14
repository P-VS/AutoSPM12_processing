function [delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles)

if ~params.func.meepi
    ppparams.subfuncstring = [ppparams.substring '_task-' ppparams.task '_bold'];
    ppparams.funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring '.json']);
    params.func.echoes = [1];
else
    ppparams.subfuncstring = [ppparams.substring '_task-' ppparams.task '_bold_e'];
    ppparams.funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring '1.json']);
end

ppparams.echoes = params.func.echoes;
ppparams.meepi = params.func.meepi;
ppparams.reorient = params.reorient;

jsondat = fileread(ppparams.funcjsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;

numdummy = floor(params.func.dummytime/tr);

if ppparams.testscript
    readvols = 20;  %Only for a quick test of the batch script
else
    readvols = Inf;
end

for ie=ppparams.echoes

    %% Load func data without dummy scans

    if ppparams.meepi
        nsubfuncstring = [ppparams.subfuncstring num2str(ie)];
    else
        nsubfuncstring = ppparams.subfuncstring;
    end

    [Vfunc,funcdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(nsubfuncstring,ppparams.subfmridir,numdummy,ppparams,readvols);

    fprintf('Write data without dummys\n')
    
    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(ppparams.funcfile{ie}, 'prefix','e');
        Vfunc(k).descrip = 'my_spmbatch - remove dummys';
        Vfunc(k).n = [k 1];
    end

    Vfunc = myspm_write_vol_4d(Vfunc,funcdat);

    ppparams.funcfile{ie} = spm_file(ppparams.funcfile{ie}, 'prefix','e');
    delfiles{numel(delfiles)+1} = {ppparams.funcfile{ie}};
    
    if ppparams.reorient   
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
    if params.func.fieldmap
        [ppparams,delfiles,keepfiles] = my_spmbatch_fieldmap(ie,ppparams,delfiles,keepfiles);
    else
        ppparams.vdm_file = '';
    end

    ppparams.prefix = '';

    %% Realignment
    if params.func.do_realignment 
        [funcdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ie,ppparams,params,keepfiles,delfiles);
        
        Vfunc = spm_vol(ppparams.funcfile{ie});
    end

    %% Topup geometric correction
    if params.func.pepolar
        [funcdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(funcdat,numdummy,ie,ppparams,params,delfiles,keepfiles);
        ppparams.prefix = ['u' ppparams.prefix];
    end
    
    %% Slice time correction
    if params.func.do_slicetime    
    
        fprintf('Do slice time correction\n')
        
        SliceTimes = jsondat.SliceTiming;
        nsl= numel(jsondat.SliceTiming);
        
        funcdat=my_spmbatch_st(funcdat,Vfunc,SliceTimes,tr);
    
        if ~params.func.meepi || contains(params.func.combination,'none')
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
        if ~params.func.meepi || contains(params.func.combination,'none')
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
    
    if params.func.meepi
        tefuncdata{ie}.data = funcdat;
        tefuncdata{ie}.Vfunc = Vfunc;
    end
end
    
%% Combine multiple TE timeseries for ME-fMRI
if params.func.meepi && ~contains(params.func.combination,'none')
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
    if params.func.do_normalization
        [wfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_normalization_func(ie,ppparams,params,delfiles,keepfiles);
    end
     
    %% Smooth func        
    if params.func.do_smoothing
        fprintf('Start smoothing \n')
    
        if ~exist('wfuncdat','var')
            [~,sfname,~] = fileparts(ppparams.funcfile{ie});

            [Vfunc,wfuncdat,ppparams.funcfile{ie}] = my_spmbatch_readSEfMRI(sfname,ppparams.subfmridir,0,ppparams,Inf);
    
        end

        Vfunc = spm_vol(ppparams.funcfile{ie});

        sfuncdat = zeros([Vfunc(1).dim(1),Vfunc(1).dim(2),Vfunc(1).dim(3),numel(Vfunc)]);

        spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');

        for i=1:numel(Vfunc)
            [pth,nm,~] = fileparts(Vfunc(i).fname);

            Q = fullfile(pth, ['s' nm  '.nii,' num2str(Vfunc(i).n)]);
            swfuncdat = my_spmbatch_smooth(wfuncdat(:,:,:,i),Vfunc(i),Q,[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

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