function roianalysis = tbx_cfg_acid_roi_analysis

% source image
maps         = cfg_files;
maps.tag     = 'maps';
maps.name    = 'Maps';
maps.help    = {'Select images/maps for ROI analysis'};
maps.filter  = 'image';
maps.ufilter = '.*';
maps.num     = [1 Inf];

% global ROIs
mask_glob         = cfg_files;
mask_glob.tag     = 'mask_glob';
mask_glob.name    = 'Global ROIs';
mask_glob.help    = {'Select binary masks in which metrics will be extracted.'};
mask_glob.filter  = 'image';
mask_glob.ufilter = '.*';
mask_glob.num     = [0 Inf];
mask_glob.val     = {{''}};

% subject-specific ROIs
mask_subj         = cfg_files;
mask_subj.tag     = 'mask_subj';
mask_subj.name    = 'Subject-specific ROIs';
mask_subj.help    = {'Select subject-specific binary masks in which metrics will be extracted.'};
mask_subj.filter  = 'image';
mask_subj.ufilter = '.*';
mask_subj.num     = [0 Inf];
mask_subj.val     = {{''}};

% reliability masks
mask_rel         = cfg_files;
mask_rel.tag     = 'mask_rel';
mask_rel.name    = 'Reliability masks';
mask_rel.help    = {'Select reliability masks.'};
mask_rel.filter  = 'image';
mask_rel.ufilter = '.*';
mask_rel.num     = [0 Inf];
mask_rel.val     = {{''}};

% dummy for slice-wise
dummy_slicewise        = cfg_menu;
dummy_slicewise.tag    = 'dummy_slicewise';
dummy_slicewise.name   = 'Slice-wise extraction';
dummy_slicewise.help   = {'Choose whether the metric extraction from the ROI'};
dummy_slicewise.labels = {'Yes','No'};
dummy_slicewise.values = {1 0};
dummy_slicewise.val    = {0};


% exbranch:
roianalysis        = cfg_exbranch;
roianalysis.tag    = 'roianalysis';
roianalysis.name   = 'ROI analysis';
roianalysis.help   = {'Extraction of average metrics from specified ROIs.'};
roianalysis.val    = {maps, mask_glob, mask_subj, mask_rel, dummy_slicewise};
roianalysis.prog   = @local_roi_analysis;

end

function out = local_roi_analysis(job)
    acid_roi_analysis(char(job.maps), char(job.mask_glob), char(job.mask_subj), char(job.mask_rel),...
        job.dummy_slicewise);
end