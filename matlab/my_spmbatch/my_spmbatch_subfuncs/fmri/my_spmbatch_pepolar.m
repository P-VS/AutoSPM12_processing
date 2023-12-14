function [outfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(infuncdat,numdummy,ne,ppparams,params,delfiles,keepfiles)

fdir = fullfile(ppparams.subpath,'fmap');

if ~ppparams.meepi
    subfmapstring = [ppparams.substring '_dir-pi_epi'];

    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(subfmapstring,fdir,numdummy,ppparams,1);

    ppfunc = fullfile(fdir,[subfmapstring '.nii']);
else%if contains(params.combination,'none')
    subfmapstring = [ppparams.substring '_dir-pi_epi_e' num2str(ne)];

    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(subfmapstring,fdir,numdummy,ppparams,1);

    ppfunc = fullfile(fdir,[subfmapstring '.nii']);
end

Vppfunc.fname = spm_file(ppfunc, 'prefix','f');
Vppfunc.descrip = 'my_spmbatch - first volume';
Vppfunc.n = [1 1];
Vppfunc = spm_create_vol(Vppfunc);
Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);

ppfunc = spm_file(ppfunc, 'prefix','f');
delfiles{numel(delfiles)+1} = {ppfunc};

if ppparams.reorient
    MMfile = fullfile(fdir,[ppparams.substring '_MM_dir-pi_epi_e' num2str(ppparams.echoes(1)) '.mat']);

    if ne==ppparams.echoes(1)
        auto_acpc_reorient(ppfunc,'EPI');

        Vppfunc = spm_vol(ppfunc);
        MM = Vppfunc(1).mat;

        save(MMfile,'MM','-mat')

        delfiles{numel(delfiles)+1} = {MMfile};
    else
        MM = load(MMfile);
        MM = MM.MM;

        Vppfunc = spm_vol(ppfunc);
        Vppfunc = my_reset_orientation(Vppfunc,MM);
    end
end

%% coregister fmap to func

CMfile = fullfile(fdir,[ppparams.substring '_CM_dir-pi_epi_e' num2str(ppparams.echoes(1)) '.mat']);

if ne==ppparams.echoes(1)
    estwrite.ref(1) = {ppparams.reffunc{ne}};
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
    
    save(CMfile,'out_coreg','-mat')

    delfiles{numel(delfiles)+1} = {CMfile};
else
    in_coreg = load(CMfile);
    in_coreg = in_coreg.out_coreg;

    write.ref = {in_coreg.rfiles{1}};
    write.source = {ppfunc};
    write.roptions.interp = 4;
    write.roptions.wrap = [0 0 0];
    write.roptions.mask = 0;
    write.roptions.prefix = 'r';

    out_coreg = spm_run_coreg(write);
end

delfiles{numel(delfiles)+1} = {ppfunc};
delfiles{numel(delfiles)+1} = {spm_file(ppfunc, 'prefix','r')};

%% HYSCO correction
    
jsondat = fileread(ppparams.funcjsonfile);
jsondat = jsondecode(jsondat);
pedir = jsondat.PhaseEncodingDirection;

if contains(pedir,'i')
    pedim = 1;
    WrapD = [1 0 0];
elseif contains(pedir,'j')
    pedim = 2;
    WrapD = [0 1 0];
else
    pedim = 3;
    WrapD = [0 0 1];
end

if contains(pedir,'-')
    blipdir=-1;
else
    blipdir=1;
end

V = spm_vol(ppparams.funcfile{ne});
for ti=1:numel(V)
    flist{ti,1} = [ppparams.funcfile{ne} ',' num2str(ti)];
end

if blipdir==1
    source_up = ppparams.reffunc{ne};
    source_dw = out_coreg.rfiles{1};

    others_up = infuncdat;
    others_dw = [];

    derfolder = fullfile(ppparams.subpath,'func');
else
    source_up = out_coreg.rfiles{1};
    source_dw = ppparams.reffunc{ne};

    others_up = [];
    others_dw = infuncdat;

    derfolder = fullfile(ppparams.subpath,'fmap');
end

d = dir(derfolder);
ddir=[d(:).isdir];
dfolders = {d(ddir).name};
delfolders = find(contains(dfolders,'derivatives'));

if ~isempty(delfolders)
    for delf=1:numel(delfolders)
        rmdir(fullfile(derfolder,dfolders{delfolders(delf)}),'s');
    end
end

dummy_fast  = 1;
alpha       = acid_get_defaults('hysco.alpha');
beta        = acid_get_defaults('hysco.beta');
restrictdim = acid_get_defaults('hysco.restrictdim');
dummy_ecc   = acid_get_defaults('hysco.dummy_ecc');
res         = acid_get_defaults('hysco.resample');

[VB1,~,~,hpth] = acid_hysco_main(source_up,source_dw,[],[],pedim,dummy_fast,dummy_ecc,alpha,beta,1,restrictdim,res,{''});

[~,nm,~] = fileparts(VB1.fname);
VB1fname = fullfile(hpth,[nm '.nii']); %correct error in HYSCO

outfuncdat = my_spmbatch_hysco_write(source_up,others_up,others_dw,VB1fname,pedim,1);

delfiles{numel(delfiles)+1} = {fullfile(derfolder,'derivatives')}; 