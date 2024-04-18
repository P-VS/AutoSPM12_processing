function resamp = tbx_cfg_acid_resample

% source images
sources         = cfg_files;
sources.tag     = 'sources';
sources.name    = 'Images to resample';
sources.help    = {'Select the dMRI dataset you want to resample. You can load in either a set of 3d nifti or a single 4d nifti file.'};
sources.filter  = 'image';
sources.ufilter = '.*';
sources.num     = [0 Inf];
    
% voxel dimensions
res         = cfg_entry;
res.tag     = 'res';
res.name    = 'Voxel size in mm';
res.help    = {'Provide a 1 x 3  - vector with the desired voxel size.'};
res.strtype = 'e';
res.num     = [1 3];
res.val     = {[1 1 1]};
    
% interpolation
interp         = cfg_entry;
interp.tag     = 'interp';
interp.name    = 'Interpolation order';
interp.help    = {'Interpolation order as defined in spm_reslice. A splin-interpolation is used by default.'};
interp.strtype = 'e';
interp.num     = [1 1];
interp.val     = {-7};

%% Exbranch
resamp       = cfg_exbranch;
resamp.tag   = 'resamp';
resamp.name  = 'Resampling';
resamp.val   = {sources res interp};
resamp.help  = {
                    'This function resamples the DTI dataset to a resolution of choice. NOTE: The header rotation matrix will be changed. The data will be resampled in dicom space.'
                    'Interpolation to a higher spatial resolution might be advantageous for improved delination of small structures in the brain.' 
                    'For spinal cord DTI (Mohammadi et al., Neuroimage, 2013), we showed that interpolation to higher in-plane spatial resolution increased the effective resolution of the tensor estimates and thus improved deliniation of the butterfly-shaped gray matter structure in the spinal cord.'}';
resamp.prog  = @local_resampletool;
resamp.vout  = @vout_resampletool;

end

function out = local_resampletool(job)

    dummy_mask = false;
    VG = acid_resample(char(job.sources), job.res, nan(2,3), dummy_mask, job.interp);    
    VG = VG(:);
    [a,b,c] = fileparts(VG(1).fname);
    
    for i = 1:size(VG,1)   
        out.rfiles{i,:} = [a filesep b c ',' num2str(i)];   
    end
end

function dep = vout_resampletool(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Resampled images';
    dep(1).src_output = substruct('.','rfiles');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end