function [outfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(infuncdat,numdummy,ne,nt,ppparams,delfiles,keepfiles)

if nt==1
    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(ppparams.subfmapdir,ppparams.func(ne).fmapfile,numdummy+1,ppparams,1);
    
    Vppfunc.fname = fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]);
    Vppfunc.descrip = 'my_spmbatch - first volume';
    Vppfunc.n = [1 1];
    Vppfunc = spm_create_vol(Vppfunc);
    Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);
    
    auto_acpc_reorient([fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]) ',1'],'EPI');
    MM = Vppfunc(1).mat;
    
    Vppfunc = my_reset_orientation(Vppfunc,MM);
    Vppfunc = spm_create_vol(Vppfunc);
    
    ppfunc = fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]);
    delfiles{numel(delfiles)+1} = {ppfunc};
    
    clear ppfuncdat Vppfunc
    
    %% coregister fmap to func
    
    reffunc = [ppparams.func(ne).tprefix ppparams.func(ne).funcfile];
    
    Tfunc = spm_vol(fullfile(ppparams.subfuncdir,[reffunc ',1']));
    Tdat = spm_read_vols(Tfunc);
    
    refsplit1 = split(reffunc,'sub-');
    refsplit2 = split(refsplit1{2},'_');
    
    nreffunc = ['sub-' refsplit2{1} '_' refsplit1{1} '_' refsplit1{2}];
    
    Tfunc.fname = fullfile(ppparams.subfuncdir,nreffunc);
    Tfunc = spm_write_vol(Tfunc,Tdat);
    
    reffunc = fullfile(ppparams.subfuncdir,nreffunc);
    
    clear Tfunc Tdat

    estwrite.ref(1) = {reffunc};
    estwrite.source(1) = {ppfunc};
    estwrite.other = {ppfunc};
    estwrite.eoptions.cost_fun = 'nmi';
    estwrite.eoptions.sep = [4 2];
    estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    estwrite.eoptions.fwhm = [7 7];
    estwrite.roptions.interp = 4;
    estwrite.roptions.wrap = [0 0 0];
    estwrite.roptions.mask = 0;
    estwrite.roptions.prefix = 'r';
    
    out_coreg = spm_run_coreg(estwrite);
    
    delfiles{numel(delfiles)+1} = {ppfunc};
    delfiles{numel(delfiles)+1} = {spm_file(ppfunc, 'prefix','r')};
    
    [pth,fname,~] = fileparts(out_coreg.rfiles{1});
    refsplit1 = split(fname,'sub-');
    refsplit2 = split(refsplit1{2},'_');
    
    ncoregf = fullfile(pth,['sub-' refsplit2{1} '_' refsplit1{1} '_' refsplit1{2} '.nii']);
    copyfile(fullfile(pth,[fname '.nii']),ncoregf);
    
    %% HYSCO correction
        
    jsondat = fileread(ppparams.func(ne).jsonfile);
    jsondat = jsondecode(jsondat);
    pedir = jsondat.PhaseEncodingDirection;
    
    if contains(pedir,'i')
        ppparams.pepolar.pedim = 1;
        WrapD = [1 0 0];
    elseif contains(pedir,'j')
        ppparams.pepolar.pedim = 2;
        WrapD = [0 1 0];
    else
        ppparams.pepolar.pedim = 3;
        WrapD = [0 0 1];
    end
    
    if contains(pedir,'-')
        ppparams.pepolar.blipdir=-1;
    else
        ppparams.pepolar.blipdir=1;
    end
    
    if ppparams.pepolar.blipdir==1
        ppparams.pepolar.source_up = reffunc;
        ppparams.pepolar.source_dw = ncoregf;
    else
        ppparams.pepolar.source_up = ncoregf;
        ppparams.pepolar.source_dw = reffunc;
    end
    
    dummy_fast  = 1;
    alpha       = acid_get_defaults('hysco.alpha');
    beta        = acid_get_defaults('hysco.beta');
    restrictdim = acid_get_defaults('hysco.restrictdim');
    dummy_ecc   = acid_get_defaults('hysco.dummy_ecc');
    res         = acid_get_defaults('hysco.resample');
    
    [ppparams.pepolar.VB1,~,~] = acid_hysco(ppparams.pepolar.source_up,ppparams.pepolar.source_dw,[],[],ppparams.pepolar.pedim,dummy_fast,dummy_ecc,alpha,beta,restrictdim,res);

    delfiles{numel(delfiles)+1} = {reffunc};
    delfiles{numel(delfiles)+1} = {ncoregf};
end

dim = size(infuncdat);
if numel(dim)<4, reshape(infuncdat,[dim(1),dim(2),dim(3),1]); end
if ppparams.pepolar.blipdir==1
    others_up = infuncdat;
    others_dw = [];
else
    others_up = [];
    others_dw = infuncdat;
end

outfuncdat = my_spmbatch_hysco_write(ppparams.pepolar.source_up,others_up,others_dw,ppparams.pepolar.VB1.dat.fname,ppparams.pepolar.pedim,1);
if numel(dim)<4, outfuncdat = reshape(outfuncdat,[dim(1),dim(2),dim(3)]); end

clear infuncdat