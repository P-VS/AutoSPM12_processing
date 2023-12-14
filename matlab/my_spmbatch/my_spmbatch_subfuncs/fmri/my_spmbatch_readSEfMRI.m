function [Vfunc,funcdat,funcfile] = my_spmbatch_readSEfMRI(subfuncstring,subfmridir,numdummy,ppparams,readvols)

%% Loading the fMRI time series and deleting dummy scans
fprintf('Reading the data \n')

funcfile = fullfile(subfmridir,[subfuncstring '.nii']);

Vfunc = spm_vol(funcfile);

if ppparams.reorient
    if ppparams.meepi
        nfname = split(subfuncstring,'bold_e');
        subfuncstring = [nfname{1} '_bold_e1'];
    end

    transfile = fullfile(subfmridir,[subfuncstring '_reorient.mat']);
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
    if readvols>1
        Vfunc = Vfunc(1:readvols);
    else
        Vfunc = Vfunc(1);
    end
end

funcdat = spm_read_vols(Vfunc);