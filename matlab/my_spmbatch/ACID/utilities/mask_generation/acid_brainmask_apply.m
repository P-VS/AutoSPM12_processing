function acid_brainmask_apply(P_mask, P, p_out)

% =========================================================================
% The script
%
%
% Created by S.Mohammadi 11/10/2012
% Adapted for multiple images by G.David
% =========================================================================
    
% load in headers for mask and image(s) to be masked
V_mask = spm_vol(P_mask);
V = spm_vol(P);
    
for i = 1:size(P,1)
        
    % load in mask and image(s) to be masked
    I_mask = acid_read_vols(V_mask, V(i), 1);  
    I = acid_read_vols(V(i), V(i), 1);
        
    % mask images
    I_masked = I.*I_mask;
        
    % define filename for masked output
    keyword = 'BRAIN-MASKED';
    fname = acid_bids_filename(V(i), keyword, 'same', '.nii');
    fname = [p_out filesep fname];
    
    % create nifti object for output
    V(i).fname = fname;
    spm_write_vol(V(i), I_masked);
    
end

end