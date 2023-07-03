function [ppparams,delfiles,keepfiles] = my_spmbatch_pepolar(numdummy,ne,ppparams,params,delfiles,keepfiles)

fdir = fullfile(ppparams.subpath,'fmap');

if ~params.meepi
    subfmapstring = [ppparams.substring '_dir-pi_epi'];

    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(subfmapstring,fdir,numdummy,params,1);

    ppfunc = fullfile(fdir,[subfmapstring '.nii']);
else
    subfmapstring = [ppparams.substring '_dir-pi_epi_e' num2str(ne)];

    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(subfmapstring,fdir,numdummy,params,1);

    ppfunc = fullfile(fdir,[subfmapstring '.nii']);
end

Vppfunc.fname = spm_file(ppfunc, 'prefix','f');
Vppfunc.descrip = 'my_spmbatch - first volume';
Vppfunc.n = [1 1];
Vppfunc = spm_create_vol(Vppfunc);
Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);

ppfunc = spm_file(ppfunc, 'prefix','f');
delfiles{numel(delfiles)+1} = {ppfunc};

if params.reorient
    MMfile = fullfile(fdir,[ppparams.substring '_MM_dir-pi_epi_e' num2str(params.echoes(1)) '.mat']);

    if ne==params.echoes(1)
        auto_acpc_reorient(ppfunc,'EPI');

        Vppfunc = spm_vol(ppfunc);
        MM = Vppfunc(1).mat;

        save(MMfile,'MM','-mat')

        delfiles{numel(delfiles)+1} = {MMfile};
    else
        MM = load(MMfile);
        MM = MM.MM;

        Vppfunc = spm_vol(ppfunc);
        Vfunc = my_reset_orientation(Vppfunc,MM);
    end
end

%% coregister fmap to func

CMfile = fullfile(fdir,[ppparams.substring '_CM_dir-pi_epi_e' num2str(params.echoes(1)) '.mat']);

if ne==params.echoes(1)
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

%% HYSCO fieldmap
    
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

hyscostep = 1;

if blipdir==1
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_up(1) = {ppparams.reffunc{ne}};
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_dw(1) = {out_coreg.rfiles{1}};
else
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_up(1) = {out_coreg.rfiles{1}};
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_dw(1) = {ppparams.reffunc{ne}};
end
matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_up(1) = {''};
matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_dw(1) = {''};
matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.perm_dim = pedim;
matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.dummy_fast = 1;
matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.outdir = {''};

delfiles{numel(delfiles)+1} = {spm_file(Vppfunc(1).fname, 'prefix','u2r')};
delfiles{numel(delfiles)+1} = {spm_file(Vppfunc(1).fname, 'prefix','HySCov2_r')};
delfiles{numel(delfiles)+1} = {spm_file(ppparams.reffunc{ne}, 'prefix','u2')};

%% Convert fieldmap into vdm

te = jsondat.EchoTime*1000;
trt = jsondat.TotalReadoutTime*1000;

fmstep = 2;

matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = cfg_dep('HySCO: Inhomogeneity field', substruct('.','val', '{}',{hyscostep}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fieldmap'));
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap = {''};
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [te te];
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdir;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = trt;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {ppparams.reffunc{ne}};
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{2}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

d = dir(fullfile(ppparams.subpath,'func'));
ddir=[d(:).isdir];
dfolders = {d(ddir).name};
delfolders = find(contains(dfolders,'derivatives'));

if ~isempty(delfolders)
    for delf=1:numel(delfolders)
        rmdir(fullfile(ppparams.subpath,'func',dfolders{delfolders(delf)}),'s');
    end
end

%%Run matlabbatch
if exist("matlabbatch",'var')
    %spm_jobman('run', matlabbatch);
end

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subpath,'func','derivatives')}; 
delfiles{numel(delfiles)+1} = {spm_file(ppparams.reffunc{ne}, 'prefix','u')};    

[pth,name,ext] = fileparts(ppparams.reffunc{ne});
ppparams.vdm_file = fullfile(ppparams.subpath,'func','derivatives','HySCO-Run',['vdm5_' name '_desc-H-HySCO-ESTIMATED-FIELDMAP_dwi.nii']);