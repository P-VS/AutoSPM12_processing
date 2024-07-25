function [rfuncdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ne,nt,nvols,ppparams,params,keepfiles,delfiles)

Vref = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile ',1']));
funcref  = spm_read_vols(Vref);

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile]));
funcdat = spm_read_vols(Vfunc(nt:nt+nvols-1));
dim = size(funcdat);
if numel(dim)<4, funcdat = reshape(funcdat,[dim(1),dim(2),dim(3),1]); end

rfuncdat = zeros([dim(1),dim(2),dim(3),nvols]);

for iv=1:nvols
    fprintf(['\nRealign vol ' num2str(nt+iv-1)])
    
    R(1,1).mat = Vref.mat;
    R(1,1).dim = Vref.dim;
    R(1,1).Vol = funcref;

    R(2,1).mat = Vfunc(iv).mat;
    R(2,1).dim = Vfunc(iv).dim;
    R(2,1).Vol = reshape(funcdat(:,:,:,iv),[Vfunc(iv).dim(1),Vfunc(iv).dim(2),Vfunc(iv).dim(3)]);

    %% estimate the realignment parameters
    if ne==params.func.echoes(1)
        tR = R;
    
        max_fr = max(tR(1,1).Vol,[],'all');
        max_fd = max(tR(2,1).Vol,[],'all');
        tR(2,1).Vol = tR(2,1).Vol * max_fr/max_fd;
    
        eoptions.quality = 0.9;
        eoptions.sep = 4;
        eoptions.fwhm = 5;
        eoptions.rtm = 0;
        eoptions.interp = 2;
        eoptions.wrap = [0 0 0];
        eoptions.PW = '';
        eoptions.lkp = 1:6;
    
        if (nt+iv-1)==1
            ppparams.realign.A0=[];
            ppparams.realign.x1=[];ppparams.realign.x2=[];ppparams.realign.x3=[];
            ppparams.realign.wt=[];
            ppparams.realign.deg=[];
            ppparams.realign.b=[];
        end
    
        % Using OpenNFT functionality: Realign to reference volume
        [tR, ppparams.realign.A0, ppparams.realign.x1, ppparams.realign.x2, ppparams.realign.x3, ...
            ppparams.realign.wt, ppparams.realign.deg, ppparams.realign.b, ~] = ...
                my_spmbatch_spm_realign_rt(tR, eoptions, nt, 1, ppparams.realign.A0, ppparams.realign.x1, ppparams.realign.x2, ppparams.realign.x3, ...
                                            ppparams.realign.wt, ppparams.realign.deg, ppparams.realign.b);
        % Get motion correction parameters
        tmpMCParam = spm_imatrix(tR(2,1).mat / tR(1,1).mat);
        if (nt+iv-1) == 1, ppparams.realign.offsetMCParam = tmpMCParam(1:6); end
        MP = tmpMCParam(1:6) - ppparams.realign.offsetMCParam;
    
        if ~isfield(ppparams,'rp_file')
            ppparams.rp_file = spm_file(fullfile(ppparams.subfuncdir,[ppparams.func(ne).tprefix ppparams.func(ne).funcfile]), 'prefix','rp_','ext','.txt');
        
            if params.func.meepi
                [rppath,rpname,~] = fileparts(ppparams.rp_file);
                nrpname = split(rpname,'_echo-');
                if contains(nrpname{2},'_aslbold'), ext='_aslbold'; else ext='_bold'; end
                nrp_file = fullfile(rppath,[nrpname{1} ext '.txt']);
        
                ppparams.rp_file = nrp_file;
            end
    
            save(ppparams.rp_file,'MP','-ascii');
        
            keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
        else
            save(ppparams.rp_file,'MP','-append','-ascii');
        end
    
        R(2,1).mat = tR(2,1).mat;
    
        ppparams.realign.R(nt+iv-1,1).mat = R(2,1).mat;
    else
        R(2,1).mat = ppparams.realign.R(nt+iv-1,1).mat;
    end
    
    %% Reslice the functional series
    roptions.quality = 0.9;
    roptions.sep = 4;
    roptions.fwhm = 5;
    roptions.which = 2;
    roptions.mean = 0;
    roptions.interp = 4;
    roptions.wrap = [0 0 0];
    roptions.mask = 1;
    roptions.prefix = 'r';

    % Reslice to reference image grid
    rfuncdat(:,:,:,iv) = my_spmbatch_spm_reslice_rt(R, roptions);

    clear R tR
end

clear funcdat funcref Vref Vfunc