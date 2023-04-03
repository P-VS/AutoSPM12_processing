function [out,norms] = denoise_acompcor(job)

warnstate = warning;
warning off;

spm_get_defaults;

rp_file = job.rp_file{1};
expand_regressors = job.expand_regressors;

confounds = load(rp_file);

if expand_regressors>1
    confounds = cat(2,confounds,cat(1,zeros(1,6),diff(confounds)));

    if expand_regressors==3
        confounds = cat(2,confounds,power(confounds,2));
    end
end

[pth,name,ext] = fileparts(job.csf_map);
split_name = split(name,'c3');

c1im = fullfile(pth,[split_name{1} 'c1' split_name{2} '.nii']);
c2im = fullfile(pth,[split_name{1} 'c2' split_name{2} '.nii']);
c3im = fullfile(pth,[split_name{1} 'c3' split_name{2} '.nii']);

GM = spm_vol(c1im);
WM = spm_vol(c2im);
CSF = spm_vol(c3im);

gmdat = spm_read_vols(GM);
wmdat = spm_read_vols(WM);
csfdat = spm_read_vols(CSF);

braindat = gmdat+wmdat;
braindat(braindat<0.2)=0;
braindat(braindat>0.0)=1;

csfdat(braindat>0.0)=0;
csfdat(csfdat<0.8)=0;
csfdat(csfdat>0.0)=1;

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

acc_confounds = den_fmri_acompcor(funcdat(:,:,:,:),{csfdat},job.ncomp,'confounds',confounds,'filter',[job.reptime job.filtfreq(1) job.filtfreq(2)],'PolOrder',1);

confounds = cat(2,confounds,acc_confounds);
    
rp_file = spm_file(rp_file, 'prefix','acc_','ext','.txt');

writematrix(confounds,rp_file,'Delimiter','tab');

norms = rp_file;

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

[funcdat,~] = den_fmri_cleaning(funcdat(:,:),1,[job.reptime job.filtfreq(1) job.filtfreq(2)],confounds,[],'restoremean','on');

funcdat = reshape(funcdat(:,:),s);

if numel(funcfile)==1
    for k=1:numel(Vfunc)
        Vfunc(k).fname = spm_file(funcf, 'prefix','d');
        Vfunc(k).descrip = 'denoise - acompcor';
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
