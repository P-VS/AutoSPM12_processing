function [ppparams,delfiles,keepfiles] = my_spmbacth_faslsubtraction(ppparams,params,delfiles,keepfiles)

GM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile));
WM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile));
CSF = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile));

gmim = spm_read_vols(GM);
wmim = spm_read_vols(WM);
csfim = spm_read_vols(CSF);

mask = gmim+wmim;
mask(mask<0.1) = 0;
mask(mask>0) = 1;

csfim(mask>0) = 0;

nfname = split(ppparams.perf(1).aslfile,'_asl');

Vasl=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix ppparams.perf(1).aslfile]));

tdim = numel(Vasl);
nvols = params.loadmaxvols;
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end
    minvol = max([1,ti-1]);
    maxvol = min([tdim,ti+nvols]);

    fprintf(['Subtract vols: ' num2str(ti) '-' num2str(minvol+nvols-1) '\n'])
    
    fasldata = spm_read_vols(Vasl(minvol:maxvol));
    
    voldim = size(fasldata);
    
    %fasldata = fasldata .* repmat(mask,1,1,1,voldim(4));

    if ti==1
        csfdata = sum(reshape(fasldata .* repmat(csfim,1,1,1,voldim(4)),[voldim(1)*voldim(2)*voldim(3),voldim(4)]),1)/numel(find(csfim>0));
        
        conidx = 2:2:tdim;
        labidx = 1:2:tdim;
        
        mean_csfcon = mean(csfdata(conidx));
        mean_csflab = mean(csfdata(labidx));
        
        if mean_csflab>mean_csfcon
            conidx = 1:2:tdim;
            labidx = 2:2:tdim;
        end
    end

    deltamdata = zeros([voldim(1),voldim(2),voldim(3),voldim(4)]);

    for iv=2:nvols-1
        if sum(conidx==(iv+ti-1))
            itcondat=fasldata(:,:,:,iv); 
        else
            itcondat=(fasldata(:,:,:,iv-1)+fasldata(:,:,:,iv+1))/2; 
        end
           
        if sum(labidx==(iv+ti-1)) 
            itlabdat=fasldata(:,:,:,iv); 
        else
            itlabdat=(fasldata(:,:,:,iv-1)+fasldata(:,:,:,iv+1))/2; 
        end
    
        deltamdata(:,:,:,iv) = itcondat-itlabdat;
    end

    if ti==1
        if sum(conidx==1), itcondat=fasldata(:,:,:,1); else itcondat=fasldata(:,:,:,2); end
        if sum(labidx==1), itlabdat=fasldata(:,:,:,1); else itlabdat=fasldata(:,:,:,2); end

        deltamdata(:,:,:,1) = itcondat-itlabdat;
    else
        deltamdata = deltamdata(:,:,:,2:end);
    end

    if maxvol==tdim
        if sum(conidx==maxvol), itcondat=fasldata(:,:,:,voldim(4)); else itcondat=fasldata(:,:,:,voldim(4)-1); end
        if sum(labidx==maxvol), itlabdat=fasldata(:,:,:,voldim(4)); else itlabdat=fasldata(:,:,:,voldim(4)-1); end

        deltamdata(:,:,:,voldim(4)) = itcondat-itlabdat;
    else
        deltamdata = deltamdata(:,:,:,1:end-1);
    end
        
    Vout = Vasl(1);
    rmfield(Vout,'pinfo');
    for iv=1:nvols
        Vout.fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii']);
        Vout.descrip = 'my_spmbatch - deltam';
        Vout.dt = [spm_type('float32'),spm_platform('bigend')];
        Vout.n = [ti+iv-1 1];
        Vout = spm_write_vol(Vout,deltamdata(:,:,:,iv));
    end

    clear deltamdata
end

ppparams.perf(1).deltamprefix = ppparams.perf(1).aslprefix;
ppparams.perf(1).deltamfile = [nfname{1} '_deltam.nii'];