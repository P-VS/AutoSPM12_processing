function [vdm_file,delfiles,keepfiles] = my_spmbatch_fieldmap(subpath,substring,task,reffunc,delfiles,keepfiles)

mbstep = 1;

%%Load fieldmap data per echo and coregister to func
    
e1dat = fullfile(subpath,'fmap',[substring '_fmap_echo-1_am.nii']);
e1json = fullfile(subpath,'fmap',[substring '_fmap_echo-1_am.json']);

if isfile(e1dat)
    e1phdat = fullfile(subpath,'fmap',[substring '_fmap_echo-1_ph.nii']);

    e1jsondat = fileread(e1json);
    e1jsondat = jsondecode(e1jsondat);
    te1 = e1jsondat.EchoTime*1000;

    Ve1amp = spm_vol(e1dat);
    Ve1ph  = spm_vol(e1phdat);
else
    e1dat = fullfile(subpath,'fmap',[substring '_fmap_echo-1.nii']);
    e1json = fullfile(subpath,'fmap',[substring '_fmap_echo-1.json']);

    e1jsondat = fileread(e1json);
    e1jsondat = jsondecode(e1jsondat);
    te1 = e1jsondat.EchoTime*1000;

    Ve1 = spm_vol(e1dat);
    Ve1amp = Ve1(1);
    Ve1ph  = Ve1(2);

    Ve1amp = do_save_intermediate_results(Ve1amp,spm_read_vols(Ve1amp),'am','spm - amplitude');
    Ve1ph = do_save_intermediate_results(Ve1ph,spm_read_vols(Ve1ph),'ph','spm - phase');

    delfiles{numel(delfiles)+1} = {Ve1amp(1).fname};
    delfiles{numel(delfiles)+1} = {Ve1ph(1).fname};
end

fme1step = mbstep;

matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {Ve1amp(1).fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other(1) = {Ve1ph(1).fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

mbstep = mbstep+1;

amp1file = spm_file(Ve1amp(1).fname, 'prefix','r');
ph1file = spm_file(Ve1ph(1).fname, 'prefix','r');

delfiles{numel(delfiles)+1} = {amp1file};
delfiles{numel(delfiles)+1} = {ph1file};

e2dat = fullfile(subpath,'fmap',[substring '_fmap_echo-2_am.nii']);
e2json = fullfile(subpath,'fmap',[substring '_fmap_echo-2_am.json']);

if isfile(e2dat)
    e2phdat = fullfile(subpath,'fmap',[substring '_fmap_echo-2_ph.nii']);

    e2jsondat = fileread(e2json);
    e2jsondat = jsondecode(e2jsondat);
    te2 = e2jsondat.EchoTime*1000;

    Ve2amp = spm_vol(e2dat);
    Ve2ph  = spm_vol(e2phdat);
else
    e2dat = fullfile(subpath,'fmap',[substring '_fmap_echo-2.nii']);
    e2json = fullfile(subpath,'fmap',[substring '_fmap_echo-2.json']);

    e2jsondat = fileread(e2json);
    e2jsondat = jsondecode(e2jsondat);
    te2 = e2jsondat.EchoTime*1000;

    Ve2 = spm_vol(e2dat);
    Ve2amp = Ve2(1);
    Ve2ph  = Ve2(2);

    Ve2amp = do_save_intermediate_results(Ve2amp,spm_read_vols(Ve2amp),'am','spm - amplitude');
    Ve2ph = do_save_intermediate_results(Ve2ph,spm_read_vols(Ve2ph),'ph','spm - phase');

    delfiles{numel(delfiles)+1} = {Ve2amp(1).fname};
    delfiles{numel(delfiles)+1} = {Ve2ph(1).fname};
end

fme2step = mbstep;

matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {Ve2amp(1).fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other(1) = {Ve2ph(1).fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

mbstep = mbstep+1;

amp2file = spm_file(Ve2amp(1).fname, 'prefix','r');
ph2file = spm_file(Ve2ph(1).fname, 'prefix','r');

delfiles{numel(delfiles)+1} = {amp2file};
delfiles{numel(delfiles)+1} = {ph2file};

%% Fieldmap

fmstep = mbstep;

funcjsonfile = fullfile(subpath,'func',[substring '_task-' task '_bold.json']);
    
jsondat = fileread(funcjsonfile);
jsondat = jsondecode(jsondat);
pedir = jsondat.PhaseEncodingDirection;

if contains(pedir,'-')
    blipdim = -1;
else
    blipdim = 1;
end

try
    trt = jsondat.TotalReadoutTime*1000;
catch
    trt = jsondat.AcquisitionMatrixPE*jsondat.EffectiveEchoSpacing*1000;
end

matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase(1) = {amp1file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag(1) = {ph1file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase(1) = {amp2file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag(1) = {ph2file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [te1 te2];
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdim;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = trt;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
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
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.session.epi(1) = {reffunc};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

mbstep = mbstep+1;

vdm_file = spm_file(amp1file, 'prefix','vdm5_sc');

delfiles{numel(delfiles)+1} = {spm_file(ph1file, 'prefix','m')};
delfiles{numel(delfiles)+1} = {spm_file(ph1file, 'prefix','bmask')};
delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','sc')};
delfiles{numel(delfiles)+1} = {spm_file(amp2file, 'prefix','sc')};
delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','fpm_sc')};
delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','vdm5_sc')};
delfiles{numel(delfiles)+1} = {spm_file(reffunc, 'prefix','u')};

%%Run matlabbatch
if exist("matlabbatch",'var')
    spm_jobman('run', matlabbatch);
end