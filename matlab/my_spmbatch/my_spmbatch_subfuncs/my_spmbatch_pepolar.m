function [vdm_file,delfiles,keepfiles] = my_spmbatch_pepolar(subpath,substring,task,numdummy,reffunc,params,delfiles,keepfiles)

mbstep = 1;

if params.nechoes==1
    fdir = fullfile(subpath,'fmap');
    subfmapstring = [substring '_dir-pi_epi'];

    [Vppfunc,ppfuncdat] = my_spmbatch_readSEfMRI(subfmapstring,fdir,numdummy,params,1);

    ppfunc = fullfile(fdir,[subfmapstring '.nii']);
else
    fdir = fullfile(subpath,'fmap');
    subfmapstring = [substring '_dir-pi_epi_e'];

    [Vppfunc,ppfuncdat,ppfunc] = my_spmbatch_readMEfMRI(subfmapstring,fdir,numdummy,params,1);
end

Vppfunc.fname = spm_file(ppfunc, 'prefix','f');
Vppfunc.descrip = 'my_spmbatch - first volume';
Vppfunc.n = [1 1];
Vppfunc = spm_create_vol(Vppfunc);
Vppfunc = spm_write_vol(Vppfunc,ppfuncdat);

ppfunc = spm_file(ppfunc, 'prefix','f');
delfiles{numel(delfiles)+1} = {ppfunc};

if params.reorient
    auto_acpc_reorient(ppfunc,'EPI');
end

%% coregister fmap to func

rfuncpsstep = mbstep;

matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {ppfunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other = {ppfunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

mbstep = mbstep+1;

delfiles{numel(delfiles)+1} = {ppfunc};
delfiles{numel(delfiles)+1} = {spm_file(ppfunc, 'prefix','r')};

%% HYSCO fieldmap
if params.nechoes==1
    funcjsonfile = fullfile(subpath,'func',[substring '_task-' task '_bold.json']);
else
    funcjsonfile = fullfile(subpath,'func',[substring '_task-' task '_bold_e1.json']);
end
    
jsondat = fileread(funcjsonfile);
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

hyscostep = mbstep;

matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_up(1) = {reffunc};
matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_dw(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{rfuncpsstep}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_up(1) = {''};
matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_dw(1) = {''};
matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.perm_dim = pedim;
matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.dummy_fast = 1;
matlabbatch{mbstep}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.outdir = {''};

mbstep = mbstep+1;

delfiles{numel(delfiles)+1} = {spm_file(Vppfunc(1).fname, 'prefix','u2r')};
delfiles{numel(delfiles)+1} = {spm_file(Vppfunc(1).fname, 'prefix','HySCov2_r')};
delfiles{numel(delfiles)+1} = {spm_file(reffunc, 'prefix','u2')};

%% Convert fieldmap into vdm

te = jsondat.EchoTime*1000;
trt = jsondat.TotalReadoutTime*1000;

if contains(pedir,'-')
    blipdir=-1;
else
    blipdir=1;
end

fmstep = mbstep;

matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = cfg_dep('HySCO: Inhomogeneity field', substruct('.','val', '{}',{hyscostep}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fieldmap'));
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap = {''};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [te te];
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdir;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = trt;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {reffunc};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

mbstep = mbstep+1;

d = dir(fullfile(subpath,'func'));
ddir=[d(:).isdir];
dfolders = {d(ddir).name};
delfolders = find(contains(dfolders,'derivatives'));

if ~isempty(delfolders)
    for delf=1:numel(delfolders)
        rmdir(fullfile(subpath,'func',dfolders{delfolders(delf)}),'s');
    end
end

[pth,name,ext] = fileparts(reffunc);
vdm_file = fullfile(subpath,'func','derivatives','HySCO-Run',['vdm5_' name '_desc-H-HySCO-ESTIMATED-FIELDMAP_dwi.nii']);

delfiles{numel(delfiles)+1} = {fullfile(subpath,'func','derivatives')}; 
delfiles{numel(delfiles)+1} = {spm_file(reffunc, 'prefix','u')};

%%Run matlabbatch
if exist("matlabbatch",'var')
    spm_jobman('run', matlabbatch);
end