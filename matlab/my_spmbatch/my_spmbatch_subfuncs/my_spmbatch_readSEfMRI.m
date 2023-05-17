function [Vfunc,funcdat] = my_spmbatch_readSEfMRI(subfuncstring,subfmridir,numdummy,params,readvols)

%% Loading the fMRI time series and deleting dummy scans
fprintf('Reading the data \n')

funcfile = fullfile(dir,[subfuncstring '.nii']);

Vfunc = spm_vol(funcfile);

if params.reorient
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
    Vfunc = Vfunc(1:readvols); %Only for a quick test of the batch script
end

funcdat = spm_read_vols(Vfunc);