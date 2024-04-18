function [outfasldata,Vout,ppparams,delfiles] = my_spmbacth_faslsubtraction(fasldata,Vasl,ppparams,params,delfiles)

GM = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(1).c1m0scanfile));
WM = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(1).c2m0scanfile));

gmim = spm_read_vols(GM);
wmim = spm_read_vols(WM);

mask = gmim+wmim;
mask_ind = find(mask>0.1);

if contains(params.asl.tagorder,'labeled')
    cidx=[-3 -2 -1 0 1 2];
    lidx=[-2 -1 0 1 2 3];
    conidx = 2:2:numel(Vasl);
    labidx = 1:2:numel(Vasl);
else
    cidx=[-2 -1 0 1 2 3];
    lidx=[-3 -2 -1 0 1 2];
    conidx = 1:2:numel(Vasl);
    labidx = 2:2:numel(Vasl);
end

voldim = size(fasldata);

outfasldata = zeros(voldim(1),voldim(2),voldim(3),voldim(4));

for iv=1:voldim(4) % sinc interpolation to replace the labeled scans by a control scan
    if isempty(find(conidx==iv)), timeshift = 2.5; else timeshift = 3; end
    idx=ceil(iv/2)+cidx;
    idx(find(idx<1))=1;
    idx(find(idx>numel(conidx)))=numel(conidx);
    nimg = fasldata(:,:,:,conidx(idx));
    nimg=reshape(nimg,size(nimg,1)*size(nimg,2)*size(nimg,3),size(nimg,4));
    clear tmpimg;
    [pn,tn]=size(nimg);
    tmpimg=sinc_interpVec(nimg(mask_ind,:),timeshift);
    Vconimg=zeros(size(nimg,1),1);
    Vconimg(mask_ind)=tmpimg;
    Vconimg=reshape(Vconimg,voldim(1),voldim(2),voldim(3));
    clear tmpimg pn tn;

    if isempty(find(labidx==iv)), timeshift = 2.5; else timeshift = 3; end
    idx=ceil(iv/2)+lidx;
    idx(find(idx<1))=1;
    idx(find(idx>numel(labidx)))=numel(labidx);
    nimg = fasldata(:,:,:,labidx(idx));
    nimg=reshape(nimg,size(nimg,1)*size(nimg,2)*size(nimg,3),size(nimg,4));
    clear tmpimg;
    [pn,tn]=size(nimg);
    tmpimg=sinc_interpVec(nimg(mask_ind,:),timeshift);
    Vlabimg=zeros(size(nimg,1),1);
    Vlabimg(mask_ind)=tmpimg;
    Vlabimg=reshape(Vlabimg,voldim(1),voldim(2),voldim(3));
    clear tmpimg pn tn;
    
    outfasldata(:,:,:,iv) = Vconimg-Vlabimg;
end

nfname = split(ppparams.asl(params.func.echoes(1)).aslfile,'_asl');

Vout = Vasl;

for j=1:numel(Vout)
    Vout(j).fname = fullfile(ppparams.subperfdir,[ppparams.asl(1).aslprefix nfname{1} '_deltam.nii']);
    Vout(j).descrip = 'my_spmbatch - deltam';
    Vout(j).pinfo = [1,0,0];
    Vout(j).dt = [spm_type('float32'),spm_platform('bigend')];
    Vout(j).n = [j 1];
end

Vout = myspm_write_vol_4d(Vout,outfasldata);

delfiles{numel(delfiles)+1} = {Vout(1).fname};

ppparams.asl(1).deltamfile = [ppparams.asl(1).aslprefix nfname{1} '_deltam.nii'];