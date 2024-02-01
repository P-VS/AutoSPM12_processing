function [delfiles,keepfiles] = my_spmbatch_cat12vbm(sub,ses,datpath,params)

ppparams.sub = sub;
ppparams.ses = ses;

ppparams.substring = ['sub-' num2str(sub,'%02d')];

if ~isfolder(fullfile(datpath,ppparams.substring))
    ppparams.substring = ['sub-' num2str(sub,'%03d')];
end

ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%03d')]);

if ~isfolder(ppparams.subpath)
    ppparams.subpath = fullfile(datpath,ppparams.substring,['ses-' num2str(ses,'%02d')]);
end

ppparams.subanatdir = fullfile(ppparams.subpath,'anat');
ppparams.preproc_anat = fullfile(ppparams.subpath,params.save_folder);

if ~exist(fullfile(ppparams.subanatdir,[ppparams.substring '_T1w_Crop_1.nii']),'file')
    nsubannat = [ppparams.substring '_T1w.nii'];
else
    nsubannat = [ppparams.substring '_T1w_Crop_1.nii'];
end
ppparams.subanat = fullfile(ppparams.subanatdir,nsubannat);

delfiles = {};
keepfiles = {};

copysubanat = spm_file(ppparams.subanat, 'prefix','r');
copyfile(ppparams.subanat,copysubanat);

[fpath,fname,~] = fileparts(copysubanat);

delfiles{numel(delfiles)+1} = {fullfile(fpath,['cat_' fname '.xml'])};
delfiles{numel(delfiles)+1} = {fullfile(fpath,['catlog_' fname '.txt'])};

catestwrite.data = {copysubanat};
catestwrite.data_wmh = {''};
catestwrite.nproc = 4;
catestwrite.useprior = '';
catestwrite.opts.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
catestwrite.opts.affreg = 'mni';
catestwrite.opts.biasacc = 0.5;
catestwrite.extopts.restypes.optimal = [1 0.3];

if params.reorient
    catestwrite.extopts.setCOM = 1;
else
    catestwrite.extopts.setCOM = 0;
end
catestwrite.extopts.APP = 1070;
catestwrite.extopts.affmod = 0;
catestwrite.extopts.spm_kamap = 0;
catestwrite.extopts.LASstr = 0.5;
catestwrite.extopts.LASmyostr = 0;
catestwrite.extopts.gcutstr = 2;
catestwrite.extopts.WMHC = 2;
catestwrite.extopts.registration.shooting.shootingtpm = {fullfile(spm('Dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_0_GS.nii')};
catestwrite.extopts.registration.shooting.regstr = 0.5;
catestwrite.extopts.vox = params.vbm.normvox;
catestwrite.extopts.bb = 12;
catestwrite.extopts.SRP = 22;
catestwrite.extopts.ignoreErrors = 1;

catestwrite.output.BIDS.BIDSno = 0;

if params.vbm.do_surface
    catestwrite.output.surface = 1;
    catestwrite.output.surf_measures = 1;

    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.central.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.central.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.pial.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.pial.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.sphere.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.sphere.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.sphere.reg.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.sphere.reg.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.white.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.white.' fname '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.pbt.' fname])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.pbt.' fname])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['lh.thickness.' fname])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['rh.thickness.' fname])};
else
    catestwrite.output.surface = 0;
    catestwrite.output.surf_measures = 0;
end

if params.vbm.do_roi_atlas
    catestwrite.output.ROImenu.atlases.neuromorphometrics = 1;
    catestwrite.output.ROImenu.atlases.lpba40 = 1;
    catestwrite.output.ROImenu.atlases.cobra = 1;
    catestwrite.output.ROImenu.atlases.hammers = 0;
    catestwrite.output.ROImenu.atlases.thalamus = 1;
    catestwrite.output.ROImenu.atlases.thalamic_nuclei = 1;
    catestwrite.output.ROImenu.atlases.suit = 1;
    catestwrite.output.ROImenu.atlases.ibsr = 0;
    catestwrite.output.ROImenu.atlases.ownatlas = {''};

    delfiles{numel(delfiles)+1} = {fullfile(fpath,['catROI_' fname '.xml'])};
    delfiles{numel(delfiles)+1} = {fullfile(fpath,['catROIs_' fname '.xml'])};

    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['catROI_' fname '.mat'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['catROIs_' fname '.mat'])};
else
    catestwrite.output.ROImenu.noROI = struct([]);
end

catestwrite.output.GM.mod = 1;
catestwrite.output.GM.dartel = 0;
catestwrite.output.WM.mod = 1;
catestwrite.output.WM.dartel = 0;
catestwrite.output.CSF.mod = 1;
catestwrite.output.CSF.dartel = 0;
catestwrite.output.CSF.warped = 0;

keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['wm' fname '.nii'])};
keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['mwp1' fname '.nii'])};
keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['mwp2' fname '.nii'])};
keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['mwp3' fname '.nii'])};

if params.vbm.do_normalization
    catestwrite.output.GM.native = 0;
    catestwrite.output.WM.native = 0;
    catestwrite.output.CSF.native = 0;
    catestwrite.output.labelnative = 0;

    delfiles{numel(delfiles)+1} = {copysubanat};

else
    catestwrite.output.GM.native = 1;
    catestwrite.output.WM.native = 1;
    catestwrite.output.CSF.native = 1;
    catestwrite.output.labelnative = 1;

    keepfiles{numel(keepfiles)+1} = {copysubanat};

    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['p0' fname '.nii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['p1' fname '.nii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['p2' fname '.nii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['p3' fname '.nii'])};
end

catestwrite.output.ct.native = 0;
catestwrite.output.ct.warped = 0;
catestwrite.output.ct.dartel = 0;
catestwrite.output.pp.native = 0;
catestwrite.output.pp.warped = 0;
catestwrite.output.pp.dartel = 0;
catestwrite.output.WMH.native = 0;
catestwrite.output.WMH.warped = 0;
catestwrite.output.WMH.mod = 0;
catestwrite.output.WMH.dartel = 0;
catestwrite.output.SL.native = 0;
catestwrite.output.SL.warped = 0;
catestwrite.output.SL.mod = 0;
catestwrite.output.SL.dartel = 0;
catestwrite.output.TPMC.native = 0;
catestwrite.output.TPMC.warped = 0;
catestwrite.output.TPMC.mod = 0;
catestwrite.output.TPMC.dartel = 0;
catestwrite.output.atlas.native = 0;
catestwrite.output.label.native = 1;
catestwrite.output.label.warped = 0;
catestwrite.output.label.dartel = 0;
catestwrite.output.bias.warped = 1;
catestwrite.output.las.native = 0;
catestwrite.output.las.warped = 0;
catestwrite.output.las.dartel = 0;
catestwrite.output.jacobianwarped = 0;
catestwrite.output.warps = [0 0];
catestwrite.output.rmat = 0;

cat_run(catestwrite);

delfiles{numel(delfiles)+1} = {fullfile(fpath,'mri')};
delfiles{numel(delfiles)+1} = {fullfile(fpath,['catlog_' fname '.txt'])};

keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['cat_' fname '.mat'])};
keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['catreport_' fname '.pdf'])};
keepfiles{numel(keepfiles)+1} = {fullfile(fpath,['catreportj_' fname '.jpg'])};