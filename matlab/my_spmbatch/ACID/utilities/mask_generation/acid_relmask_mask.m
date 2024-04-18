function acid_relmask_mask(rmse, thr)

% ================================================================================
% The function creates binary reliability mask (M) by comparing the model-fit error (e) 
% with the threshold value (thr) in each voxel. Voxels with model-fit error above the
% threshold will be masked out:
%
%   M(x) = 1, if e(x)<=thr
%   M(x) = 0, if e(x)>thr
%
% Inputs:
%   rmse  - model-fit error maps of a single or multiple subjects;
%                the script expects a char array with the filenames of
%                model-fit error maps in the rows;
%                use spm_select to generate this input
%
%   thr   - model-fit error threshold value for reliability masking;
%           this value can be either set to an arbitrary value or an optimal
%           value can be obtained using the script '';
%           no default value
%
%   p_out - output folder, default: folder of input images
%
% Outputs:
%   - binary mask for each subject specified in the input; these masks will
%     be saved in the folder of the input files
%
% Created by G. David, May 2017
% ================================================================================

% check input
if numel(rmse)==0
    error('No model-fit error maps specified.');
elseif numel(thr)==0
    error('No threshold value specified.');
end
if thr<=0
    error('The specified threshold has to be at least 0.');
end

keyword = 'REL-MASK';

for k = 1:size(rmse,1)
    
    % load in model-fit error map
    V_rmse = spm_vol(rmse(k,:));
    err = acid_read_vols(V_rmse,V_rmse,1);
    
    fprintf('Processing subject: %s\n', V_rmse.fname);
        
    % create binary reliability mask
    mask_rel = double(err <= thr);
    
    % define output directory

        [path,~,~] = spm_fileparts(V_rmse.fname);
        pout = path;

    
    % generate filename
    fname = acid_bids_filename(V_rmse, keyword, '', '.nii');
    fname = [pout filesep fname];
    
    % save reliability mask
    V_rmse.fname = fname;
    spm_write_vol(V_rmse, mask_rel);
    
end

disp('Job done!')

end