function relmask = tbx_cfg_acid_relmask_mask

% model-fit error maps
rmse         = cfg_files;
rmse.tag     = 'rmse';
rmse.name    = 'Model-fit error maps';
rmse.help    = {'Select the model-fit error maps.'};
rmse.filter  = 'image';
rmse.ufilter = '.*';
rmse.num     = [0 Inf];

% threshold for reliability masking
thr         = cfg_entry;
thr.tag     = 'thr';
thr.name    = 'Threshold for reliability masking';
thr.help    = {'The threshold value in model-fit error above which voxels are rejected.'
                   'This can be either set to an arbitrary value or to a value determined automatically.'};
thr.strtype = 'e';
thr.num     = [1 1];
thr.val     = {};

% exbranch
relmask      = cfg_exbranch;
relmask.tag  = 'relmask';
relmask.name = 'Create reliability masks';
relmask.help = {'Post-processing'};
relmask.val  = {rmse, thr};
relmask.prog = @local_REL_create_mask;
relmask.vout = @vout_REL_create_mask;

end

function out = local_REL_create_mask(job)
    acid_relmask_mask(char(job.RMS_error), job.thr);
    
%    out.masks = 
    out.masks = acid_select(out.masks(:));
    out.masks = out.masks(:);
end

function dep = vout_REL_create_mask(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Reliability masks';
    dep(1).src_output = substruct('.','masks');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end