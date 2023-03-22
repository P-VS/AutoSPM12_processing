function [out] = denoise_bppassfilt(job)

warnstate = warning;
warning off;

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

[funcdat,~] = den_fmri_cleaning(funcdat(:,:),1,[job.reptime job.filtfreq(1) job.filtfreq(2)],[],[],'restoremean','on');

funcdat = reshape(funcdat(:,:),s);

if numel(funcfile)==1
    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(funcf, 'prefix','d');
        Vfunc(k).descrip = 'denoise - bandpass filter';
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
        Vfunc.descrip = 'denoise - bandpass filter';
        Vfunc = spm_create_vol(Vfunc(k));
        Vfunc = spm_write_vol(Vfunc,funcdat(:,:,:,k));
        out = [out;spm_file(funcfile{k}, 'prefix','d')];
    end
end

end