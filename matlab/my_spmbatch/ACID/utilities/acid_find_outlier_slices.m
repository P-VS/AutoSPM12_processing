function [idx] = acid_find_outlier_slices(P,s)

%==========================================================================
% The scripts is designed to find outlier slices within a dMRI dataset.
%
% Inputs:
%   P       - dMRI dataset
%   s       - scaling factor determining the sensitivity of the outlier
%             detection
% Outpust:
%   idx     - list of outlier slices in 2d form (first row: index of volume, second row: index of slice)
%
% Created by G.David
%==========================================================================

% check input
if ~exist('P','var')
    P = char(cfg_getfile(Inf,'IMAGE','Select images','','','.*'));
end

if ~exist('s','var')
    s = 1.5;
end

% initialization
V   = spm_vol(P(1,:));
tmp = spm_read_vols(V); 
I   = zeros([size(tmp),size(P,1)]);
REF = zeros(size(tmp,3),size(P,1));

% load in the whole dMRI dataset
V = spm_vol(P);
for k = 1:size(P,1)
    I(:,:,:,k) = spm_read_vols(V(k));   
end

% create a difference matrix for each slie
I0 = squeeze(mean(mean(I,1),2));
I1 = permute(circshift(permute(I0,[2 1]),1),[2 1]);  I1(:,1) = I0(:,1);
I2 = permute(circshift(permute(I0,[2 1]),-1),[2 1]); I2(:,end) = I0(:,end);
REF(:,1)       = I2(:,1);
REF(:,2:end-1) = (I1(:,2:end-1) + I2(:,2:end-1))/2;
REF(:,end)     = I1(:,end);
DIFF = abs(I0-REF);

% remove first and last column
DIFF(:,[1,end]) = [];

% determine threshold
thr = s * std(DIFF(:));

% extract slices
[i,j] = ind2sub(size(DIFF),find(DIFF>=thr));
idx(1,:) = j'+1; % index of volumes
idx(2,:) = i'; % index of slices

% save indices
save('outliers.mat','ind')

end