function [outfuncdat,ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(infuncdat,numdummy,ne,ppparams,params,delfiles,keepfiles)

[Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(ppparams.subfmapdir,ppparams.func(ne).fmapfile,numdummy,ppparams,1);

Vppfunc.fname = fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]);
Vppfunc.descrip = 'my_spmbatch - first volume';
Vppfunc.n = [1 1];
Vppfunc = spm_create_vol(Vppfunc);
Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);

ppfunc = fullfile(ppparams.subfmapdir,['f' ppparams.func(ne).fmapfile]);
delfiles{numel(delfiles)+1} = {ppfunc};

%% coregister fmap to func

if isfield(ppparams.func(ne),'reffuncfile')
    reffunc = fullfile(ppparams.subfuncdir,ppparams.func(ne).reffuncfile);
else
    reffunc = fullfile(ppparams.subfuncdir,[ppparams.func(ne).funcfile ',1']);
end

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

%% HYSCO correction
    
jsondat = fileread(ppparams.func(ne).jsonfile);
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

V = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).prefix ppparams.func(ne).funcfile]));
for ti=1:numel(V)
    flist{ti,1} = [V(1).fname ',' num2str(ti)];
end

if blipdir==1
    source_up = reffunc;
    source_dw = out_coreg.rfiles{1};

    others_up = infuncdat;
    others_dw = [];

    derfolder = fullfile(ppparams.subpath,'func');
else
    source_up = out_coreg.rfiles{1};
    source_dw = reffunc;

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