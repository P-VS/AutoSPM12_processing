function [Vfunc,funcdat,funcfile] = my_spmbatch_readSEfMRI(subfuncstring,subfmridir,numdummy,params,readvols)

%% Loading the fMRI time series and deleting dummy scans
fprintf('Reading the data \n')

funcfile = fullfile(subfmridir,[subfuncstring '.nii']);

Vfunc = spm_vol(funcfile);

if params.reorient
    if params.meepi
        nfname = split(subfuncstring,'bold_e');
        subfuncstring = [nfname{1} '_bold_e1'];
    end

    transfile = fullfile(subfmridir,[subfuncstring '_reorient.mat']);
    if isfile(transfile)
        load(transfile,'M')
        transM = M;
    else
        transM = eye(4);
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