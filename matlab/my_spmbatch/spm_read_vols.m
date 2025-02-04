function [Y,XYZmm] = spm_read_vols(V,mask)
% Read in entire image volumes
% FORMAT [Y,XYZmm] = spm_read_vols(V,mask)
% V      - vector of mapped image volumes to read in (from spm_vol)
% mask   - implicit zero mask?
%
% Y      - 4D matrix of image data, fourth dimension indexes images
% XYZmm  - 3xn matrix of XYZ locations returned {mm}
%__________________________________________________________________________
%
% For image data types without a representation of NaN (see spm_type),
% implicit zero masking can be used. If mask is set, then zeros are
% treated as masked, and returned as NaN.
%__________________________________________________________________________
% Copyright (C) 1999-2016 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_read_vols.m 6918 2016-11-02 14:33:11Z guillaume $


%-Argument checks
%--------------------------------------------------------------------------
if nargin<2, mask = false; end

spm_check_orientations(V);

%-Read in image data
%--------------------------------------------------------------------------
n = numel(V);                                           %-#images
Y = zeros([V(1).dim(1:3),n]);                           %-image data matrix

Range = [1:n];
for i=1:n
%    for p=1:V(1).dim(3)
%        Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
%    end
    ni=V(i).n;
    Range(i) = ni(1);
end

Nii = my_spmbatch_load_nii(V(1).fname,Range);
tY = Nii.img;

for i=1:n
    ttY = tY(:,:,:,i);
    ttY = ttY * V(i).pinfo(1,1) + V(i).pinfo(2,1);
    
    Y(:,:,:,i) = ttY;
end

%-Apply implicit zero mask for image datatypes without a NaNrep
%--------------------------------------------------------------------------
if mask
    %-Work out images without NaNrep
    im = false(n,1);
    for i=1:n, im(i)=~spm_type(V(i).dt(1),'NaNrep'); end
    %-Mask
    Y(Y(:,:,:,im)==0) = NaN;
end

%-Return as 3D matrix if single image
%--------------------------------------------------------------------------
if n==1, Y=Y(:,:,:,1); end

%-Compute XYZmm co-ordinates (if required)
%--------------------------------------------------------------------------
if nargout > 1
    [R,C,P] = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
    RCP     = [R(:)';C(:)';P(:)';ones(1,numel(R))];
    XYZmm   = V(1).mat(1:3,:)*RCP;
end
