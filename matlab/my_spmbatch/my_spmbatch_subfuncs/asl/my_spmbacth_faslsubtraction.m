function [ppparams,delfiles,keepfiles] = my_spmbacth_faslsubtraction(ppparams,params,delfiles,keepfiles)

GM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile));
WM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile));
CSF = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile));

gmim = spm_read_vols(GM);
wmim = spm_read_vols(WM);
csfim = spm_read_vols(CSF);

csfim(gmim+wmim>0) = 0;
csfim(csfim<0.2) = 0;
csfim(csfim>0) = 1;

jsondat = fileread(ppparams.func(1).jsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;

Vasl=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix ppparams.perf(1).aslfile]));
fasldata = spm_read_vols(Vasl);

voldim = size(fasldata);

Vlabel=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).labelprefix ppparams.perf(1).labelfile]));
labeldata = spm_read_vols(Vlabel);
labeldata = labeldata - repmat(mean(labeldata,4),[1,1,1,voldim(4)]);

fasldata = fasldata + labeldata;

mask = my_spmbatch_mask(fasldata);

csfdata = sum(reshape(fasldata .* repmat(csfim,[1,1,1,voldim(4)]),[voldim(1)*voldim(2)*voldim(3),voldim(4)]),1)/numel(find(csfim>0));

conidx = 2:2:voldim(4);
labidx = 1:2:voldim(4);

mean_csfcon = mean(csfdata(conidx));
mean_csflab = mean(csfdata(labidx));

if mean_csflab>mean_csfcon
    conidx = 1:2:voldim(4);
    labidx = 2:2:voldim(4);
end

ppparams.asl.conidx = conidx;
ppparams.asl.labidx = labidx;

clear csfdata mean_csflab mean_csfcon

deltamdata = zeros([voldim(1)*voldim(2)*voldim(3),voldim(4)]);

fasldata = reshape(fasldata,[voldim(1)*voldim(2)*voldim(3),voldim(4)]);

condat = fasldata(mask>0,conidx);
labdat = fasldata(mask>0,labidx);

ncondat = spline((conidx-1)*tr,condat,[0:tr:(voldim(4)-1)*tr]);
nlabdat = spline((labidx-1)*tr,labdat,[0:tr:(voldim(4)-1)*tr]);

condat = fasldata(mask>0,conidx);
labdat = fasldata(mask>0,labidx);

deltamdata(mask>0,:) = ncondat-nlabdat;

deltamdata = reshape(deltamdata,[voldim(1),voldim(2),voldim(3),voldim(4)]);

clear ncondat nlabdat condat labdat fasldata

nfname = split(ppparams.perf(1).aslfile,'_asl');

if contains(params.asl.temp_resolution,'original')
    Vout = Vasl(1);
    rmfield(Vout,'pinfo');
    for iv=1:voldim(4)
        Vout.fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii']);
        Vout.descrip = 'my_spmbatch - deltam';
        Vout.dt = [spm_type('float32'),spm_platform('bigend')];
        Vout.n = [iv 1];
        Vout = spm_write_vol(Vout,deltamdata(:,:,:,iv));
    end
elseif contains(params.asl.temp_resolution,'reduced')
    nvols = floor(params.asl.dt/tr); 
    
    rdeltamdata = zeros([voldim(1),voldim(2),voldim(3),ceil(voldim(4)/nvols)]);

    Vout = Vasl(1);
    rmfield(Vout,'pinfo');
    for iv=1:ceil(voldim(4)/nvols)
        minvol=max(2,(iv-1)*nvols+1);
        maxvol=min(iv*nvols,voldim(4)-1);

        rdeltamdata = mean(deltamdata(:,:,:,minvol:maxvol),4);

        Vout.fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii']);
        Vout.descrip = 'my_spmbatch - deltam';
        Vout.dt = [spm_type('float32'),spm_platform('bigend')];
        Vout.n = [iv 1];
        Vout = spm_write_vol(Vout,rdeltamdata);

        clear rdeltamdata
    end
elseif contains(params.asl.temp_resolution,'only_mean')
    mdeltamdata = mean(deltamdata(:,:,:,2:voldim(4)-1),4);

    Vout = Vasl(1);
    rmfield(Vout,'pinfo');
    Vout.fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii']);
    Vout.descrip = 'my_spmbatch - deltam';
    Vout.dt = [spm_type('float32'),spm_platform('bigend')];
    Vout.n = [1 1];
    Vout = spm_write_vol(Vout,mdeltamdata);

    clear mdeltam
end

clear deltamdata

ppparams.perf(1).deltamprefix = ppparams.perf(1).aslprefix;
ppparams.perf(1).deltamfile = [nfname{1} '_deltam.nii'];

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii'])};