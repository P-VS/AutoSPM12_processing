function I = acid_read_vols(V, VG, interp, z)

% =========================================================================
% read image volume 
% FORMAT Aout = acid_read_vols(V,VG,res,z)
% 
% Input:
%   V      - structure containing image volume information of ith image
%   VG     - structure containing image volume information of reference image
%   interp - resampling function / order of sinc (-7 recommended)
%   z      - z position
%
% Output:
%   I   - 2D or 3D image in the space of the target image VG
%
% S.Mohammadi 06.05.2020
% =========================================================================

    dm = VG.dim;
    if exist('z','var')
        M = VG.mat*spm_matrix([0 0 z]);
        I = spm_slice_vol(V,V.mat\M,dm(1:2),interp);
    else
        I = zeros(dm);
        for z = 1:dm(3)
            M = VG.mat*spm_matrix([0 0 z]);
            I(:,:,z) = spm_slice_vol(V,V.mat\M,dm(1:2),interp);
        end
    end
end