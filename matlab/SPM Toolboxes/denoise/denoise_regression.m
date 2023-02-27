function [out,norms] = denoise_regression(job)

warnstate = warning;
warning off;

spm_defaults;

rp_file = job.rp_file{1};
expand_regressors = job.expand_regressors;

confounds = load(rp_file);

if expand_regressors>1
    confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));

    if expand_regressors==3
        confounds = cat(2,confounds,power(confounds,2));
    end

    rp_file = spm_file(rp_file, 'prefix','der_','ext','.txt');

    writematrix(confounds,rp_file,'Delimiter','tab');
end

norms = rp_file;

funcfile = job.func_file;

if numel(funcfile)==1
    [pth,name,ext] = fileparts(funcfile{1});
    funcf = fullfile(pth,[name '.nii']);

    Vfunc = spm_vol(funcf);
    funcdat = spm_read_vols(Vfunc);
else
    Vfunc = spm_vol(funcfile{1});
    dim = V.dim;

    funcdat = zeros([dim(1:3) numel(funcfile)]);

    for ti=1:numel(funcfile)
        V = spm_vol(funcfile{ti});
        funcdat(:,:,:,ti) = spm_read_vols(V);
    end
end

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

[funcdat,~] = den_fmri_cleaning(funcdat(:,:),2,[],confounds,[],'restoremean','on');

funcdat = reshape(funcdat(:,:),s);

if numel(funcfile)==1
    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(funcf, 'prefix','d');
        Vfunc(k).descrip = 'denoise - regression';
        Vfunc(k).n = [k 1];
        Vfunc(k) = spm_create_vol(Vfunc(k));
        Vfunc(k) = spm_write_vol(Vfunc(k),funcdat(:,:,:,k));
    end

    out = [spm_file(funcf, 'prefix','d')];
else
    out=[];
    for k=1:numel(funcfile)
        Vfunc=spm_vol(funcfile{k});
        Vfunc.fname = spm_file(funcfile{k}, 'prefix','d');
        Vfunc.descrip = 'denoise - regression';
        Vfunc = spm_create_vol(Vfunc(k));
        Vfunc = spm_write_vol(Vfunc,funcdat(:,:,:,k));
        out = [out;spm_file(funcfile{k}, 'prefix','d')];
    end
end

end
