function [Vfunc,funcdat] = my_spmbatch_readSEfMRI(directory,file,numdummy,ppparams,readvols)

%% Loading the fMRI time series and deleting dummy scans
fprintf('Reading the data \n')

funcfile = fullfile(directory,file);

Vfunc = spm_vol(funcfile);

if ppparams.reorient
    if contains(file,'_bold')
        nfname = split(ppparams.func(1).funcfile,'.nii');
    elseif contains(file,'_epi')
        nfname = split(ppparams.func(1).fmapfile,'.nii');
    end

    transfile = fullfile(directory,[nfname{1} '_reorient.mat']);
    if isfile(transfile)
        load(transfile,'M')
        transM = M;
    else
        transM = my_spmbatch_vol_set_com(funcfile);
        transM(1:3,4) = -transM(1:3,4);
    end

    MM = Vfunc(1).private.mat0;

    Vfunc = my_reset_orientation(Vfunc,transM * MM);
end

Vfunc = Vfunc(numdummy+1:end);

if readvols<Inf
    if readvols>numel(Vfunc), readvols = numel(Vfunc); end

    if readvols>1
        Vfunc = Vfunc(1:readvols);
    else
        Vfunc = Vfunc(1);
    end
end

funcdat = spm_read_vols(Vfunc);